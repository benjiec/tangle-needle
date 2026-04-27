from .sequence import extract_subsequence_strand_sensitive
from .match import Match, ProteinHit
from .detect import hmm_search_genome

VERBOSE_EXPAND = 0


def find_more_matches_at_locus(
  query_accession,
  hmm_file,
  old_start,
  old_end,
  target_full_sequence,
  target_accession,
  target_left,
  target_right,
  strand,
  cpus=None,
  hmm_rows=None):

    if hmm_rows is None:
        hmm_rows = hmm_search_genome(
            hmm_file, {target_accession: target_full_sequence},
            target_accession = target_accession,
            target_left = target_left,
            target_right = target_right,
            strand = strand,
            cpus = cpus,
            filter_by_evalue_cond = True  # we already assume there's a protein here...
        )

    if not hmm_rows:
        if VERBOSE_EXPAND:
            print("no hmmsearch results")
        return None

    new_matches = []
    for hmm_row in hmm_rows:
        target_sequence = extract_subsequence_strand_sensitive(target_full_sequence, hmm_row["ali_from"], hmm_row["ali_to"])
        match = Match(
            query_accession=query_accession,
            target_accession=target_accession,
            query_start=hmm_row["hmm_from"],
            query_end=hmm_row["hmm_to"],
            target_start=hmm_row["ali_from"],
            target_end=hmm_row["ali_to"],
            e_value=hmm_row["dom_evalue"],
            identity=None,
            target_sequence=target_sequence,
            matched_sequence=hmm_row["matched_sequence"]
        )
        assert match.matched_sequence == match.target_sequence_translated()
        new_matches.append(match)

    index_of_old_start = None
    index_of_old_end = None
    index_of_new_start = None
    index_of_new_end = None

    if strand > 0:
        new_matches = sorted(new_matches, key=lambda m: m.target_start)
        for i, match in enumerate(new_matches):
            if match.target_end > old_start and index_of_old_start is None:
                index_of_old_start = i
            if match.target_end >= old_end and index_of_old_end is None:
                index_of_old_end = i
        if index_of_old_end is None and index_of_old_start is not None:
            index_of_old_end = len(new_matches)-1

        if index_of_old_start is None or index_of_old_end is None or index_of_old_end < index_of_old_start:
            if VERBOSE_EXPAND:
                print("cannot find old start and end indices on fwd strand")
                for i, match in enumerate(new_matches):
                    mark = "   "
                    if i in (index_of_old_start, index_of_old_end):
                        mark = " ->"
                    if i in (index_of_new_start, index_of_new_end):
                        mark = " =>"
                    print(f" {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")
            return None

    else:
        new_matches = sorted(new_matches, key=lambda m: -m.target_start)
        for i, match in enumerate(new_matches):
            if match.target_end < old_start and index_of_old_start is None:
                index_of_old_start = i
            elif match.target_end <= old_end and index_of_old_end is None:
                index_of_old_end = i
        if index_of_old_end is None and index_of_old_start is not None:
            index_of_old_end = len(new_matches)-1

        if index_of_old_start is None or index_of_old_end is None or index_of_old_end < index_of_old_start:
            if VERBOSE_EXPAND:
                print("cannot find old start and end indices on rev strand")
                for i, match in enumerate(new_matches):
                    mark = "   "
                    if i in (index_of_old_start, index_of_old_end):
                        mark = " ->"
                    if i in (index_of_new_start, index_of_new_end):
                        mark = " =>"
                    print(f" {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")
            return None
  
    # sanity check - if we cannot collate what we think are the old matches,
    # that means what we think are the old matches are likely repeats of the
    # same match so we can just stop

    if not ProteinHit.can_collate_from_matches(new_matches[index_of_old_start:index_of_old_end+1]):
        if VERBOSE_EXPAND:
            print("sticking with the old matches", index_of_old_start, index_of_old_end)
        return None
 
    # move starting index as far back as we can
    index_of_new_start = index_of_old_start
    i = index_of_old_start-1
    while i >= 0:
        if ProteinHit.can_collate_from_matches(new_matches[i:index_of_old_end+1]):
            index_of_new_start = i
            i -= 1
        else:
            break
    
    # move ending index as far forward as we can
    index_of_new_end = index_of_old_end
    i = index_of_old_end+1
    while i <= len(new_matches)-1:
        if ProteinHit.can_collate_from_matches(new_matches[index_of_old_start:i+1]):
            index_of_new_end = i
            i += 1
        else:
            break

    if VERBOSE_EXPAND:
        print("found:")
        for i, match in enumerate(new_matches):
            mark = "  "
            if i in (index_of_old_start, index_of_old_end):
                mark = "->"
            if i in (index_of_new_start, index_of_new_end):
                mark = "=>"
            print(f"  {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")

    return new_matches[index_of_new_start:index_of_new_end+1]


def hmm_expand_protein(protein_hit, genomic_sequence_dict, hmm_file, cpus = None):

    target_full_sequence = genomic_sequence_dict[protein_hit.target_accession]
    query_accession = protein_hit.matches[0].query_accession
    target_accession = protein_hit.matches[0].target_accession
    strand = -1 if protein_hit.matches[0].on_reverse_strand else 1
    start = protein_hit.target_start
    end = protein_hit.target_end

    max_search_distance = 30000
    target_left = max(min(start, end) - max_search_distance, 1)
    target_right = min(max(start, end) + max_search_distance, len(target_full_sequence))

    if VERBOSE_EXPAND:
        print(f"{query_accession} on {target_accession}, {target_left}-{target_right} (based on {start}-{end}), strand {strand}, contig {len(target_full_sequence)}")
        for i, match in enumerate(protein_hit.matches):
            print(f" old {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")

    new_matches = find_more_matches_at_locus(
        query_accession, hmm_file, start, end,
        target_full_sequence, target_accession, target_left, target_right, strand,
        cpus = cpus
    )

    if new_matches is None:
        return protein_hit

    new_pm = ProteinHit(
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=min(m.target_start for m in new_matches) if protein_hit.target_start < protein_hit.target_end else max(m.target_start for m in new_matches),
        target_end=max(m.target_end for m in new_matches) if protein_hit.target_start < protein_hit.target_end else min(m.target_end for m in new_matches),
        hmm_file=hmm_file
    )

    return new_pm


def hmm_expand(protein_hits, genomic_sequence_dict, hmm_collection, cpus = None):
    new_protein_hits = {}

    skipped = []
    for pm in protein_hits:
        """
        print()
        print(pm.protein_hit_id)
        for nm in sorted(pm.matches, key=lambda m: m.target_start):
            print("    ", nm.target_start, nm.target_end, nm.query_start, nm.query_end)
        """

        hmm_profile = hmm_collection.get(pm.query_accession)
        if hmm_profile is None:
            if pm.query_accession not in skipped:
                skipped.append(pm.query_accession)
                print("skipping", pm.query_accession, "cannot find HMM profile")
        else:
            pm = hmm_expand_protein(pm, genomic_sequence_dict, hmm_profile, cpus = cpus)
            new_protein_hits[pm.protein_hit_id] = pm

    return list(new_protein_hits.values())

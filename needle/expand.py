from .sequence import extract_subsequence_strand_sensitive
from .match import Match, ProteinHit, build_graph, graph_paths
from .detect import hmm_search_genome

VERBOSE_EXPAND = 0


def print_matches(matches, index_of_old_start, index_of_old_end, index_of_new_start, index_of_new_end):
    for i, match in enumerate(matches):
        mark = "   "
        if i in (index_of_old_start, index_of_old_end):
            mark = " ->"
        if i in (index_of_new_start, index_of_new_end):
            mark = " =>"
        print(f" {mark} {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}, {match.e_value}, {match.score}")


def find_more_matches_at_locus(
  query_accession,
  hmm_file,
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
        if VERBOSE_EXPAND > 0:
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
            target_sequence=target_sequence,
            matched_sequence=hmm_row["matched_sequence"],
            score=hmm_row["dom_score"]
        )
        assert match.matched_sequence == match.target_sequence_translated()
        new_matches.append(match)

    return new_matches


def hmm_expand_protein(matches, genomic_sequence_dict, hmm_file, threshold = None, cpus = None, hmm_rows = None):

    target_full_sequence = genomic_sequence_dict[matches[0].target_accession]
    query_accession = matches[0].query_accession
    target_accession = matches[0].target_accession
    strand = -1 if matches[0].on_reverse_strand else 1

    if strand > 0:
        start = min([m.target_start for m in matches])
        end = max([m.target_end for m in matches])
    else:
        start = max([m.target_start for m in matches])
        end = min([m.target_end for m in matches])

    max_search_distance = 30000
    target_left = max(min(start, end) - max_search_distance, 1)
    target_right = min(max(start, end) + max_search_distance, len(target_full_sequence))

    if VERBOSE_EXPAND > 1:
        print(f"{query_accession} on {target_accession}, {target_left}-{target_right} (based on {start}-{end}), strand {strand}, contig {len(target_full_sequence)}")
        for i, match in enumerate(protein_hit.matches):
            print(f" old {match.target_start}, {match.target_end}, {match.query_start}, {match.query_end}")

    new_matches = find_more_matches_at_locus(
        query_accession, hmm_file, target_full_sequence, target_accession, target_left, target_right, strand,
        cpus = cpus,
        hmm_rows = hmm_rows
    )

    if new_matches is None:
        return []

    graph = build_graph(new_matches)
    paths = graph_paths(graph)
    proteins = []

    for path in paths:

        if threshold:
            match_total_score = sum([m.score for m in path])
            if match_total_score < threshold:
                if VERBOSE_EXPAND > 0:
                    print(f"{query_accession} on {target_accession}, {target_left}-{target_right}, score {match_total_score} not meeting threshold {threshold}")
                continue

        protein = ProteinHit(
            matches=path,
            query_start=min(m.query_start for m in path),
            query_end=max(m.query_end for m in path),
            target_start=min(m.target_start for m in path) if strand > 0 else max(m.target_start for m in path),
            target_end=max(m.target_end for m in path) if strand > 0 else min(m.target_end for m in path),
            hmm_file=hmm_file
        )
        proteins.append(protein)

    return proteins


def hmm_expand(clusters, genomic_sequence_dict, hmm_collection, thresholds = None, cpus = None):
    new_protein_hits = {}

    skipped = []
    for cluster in clusters:
        cluster_query_accession = cluster[0].query_accession
        hmm_profile = hmm_collection.get(cluster_query_accession)
        if hmm_profile is None:
            if cluster_query_accession not in skipped:
                skipped.append(cluster_query_accession)
                print("skipping", cluster_query_accession, "cannot find HMM profile")
        else:
            threshold = None
            if thresholds and cluster_query_accession in thresholds:
                threshold = thresholds[cluster_query_accession]
                if VERBOSE_EXPAND > 1:
                    print("using", threshold, "as threshold for", cluster_query_accession)
            proteins = hmm_expand_protein(cluster, genomic_sequence_dict, hmm_profile, threshold = threshold, cpus = cpus)
            for protein in proteins:
                new_protein_hits[protein.protein_hit_id] = protein

    return list(new_protein_hits.values())

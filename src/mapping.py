import mappy as mp
import pandas as pd
import logging


async def perform_alignment_async(
    fastx_file: str, ref_file: str
) -> pd.DataFrame | None:
    """
    Perform alignment between reads and the reference using Minimap2.

    :param fastx_file: Path to the FASTX file containing reads.
    :param ref_file: Path to the reference file.
    :return: PAF DataFrame if alignment is successful, None otherwise.
    """
    try:
        # Simulate alignment process
        paf_dict = perform_alignment(ref_file, str(fastx_file))
        return pd.DataFrame(paf_dict)
    except Exception as e:
        logging.error(f"Alignment failed: {e}")
        return None


def perform_alignment(reference_file: str, query_files: str) -> list:
    aligner = mp.Aligner(reference_file, preset="map-ont")  # load reference

    if not aligner:
        raise Exception("ERROR: failed to load/build index")
    results = []

    logging.debug(f"Performing alignment to {reference_file}")
    logging.debug(f"Performing alignment with  {query_files}")

    for name, seq, qual in mp.fastx_read(query_files):
        # print(seq)
        for hit in aligner.map(seq):  # traverse alignments
            results.append(
                {
                    "query_name": str(name),
                    "query_length": len(seq),
                    "start_query": hit.q_st,
                    "end_query": hit.q_en,
                    "strand": "+" if hit.strand == 1 else "-",
                    "target_name": hit.ctg,
                    "target_length": hit.ctg_len,
                    "start_target": hit.r_st,
                    "end_target": hit.r_en,
                    "aln_length": hit.blen,
                    "aln_mlength": hit.mlen,
                    "mapping_quality": hit.mapq,
                    "cigar": hit.cigar_str,
                }
            )
    logging.debug(f"{len(results)} alignments were found")
    return results


def write_paf_file(results, output_file) -> None:
    """
    takes the results of the alignment and writes a paf to disk if needed
    """
    with open(output_file, "w") as f:
        for result in results:
            f.write(
                f"{result['query_name']}\t{result['query_length']}\t{result['start_query']}\t{result['end_query']}\t{'+' if result['strand'] == '+' else '-'}\t"
                f"{result['target_name']}\t{result['target_length']}\t{result['start_target']}\t{result['end_target']}\t{result['mapping_quality']}\t{result['cigar']}\t\n"
            )

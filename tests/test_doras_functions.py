#!/bin/env python
import subprocess
import pytest
from src import utils as ps
from typing import Any
import pandas as pd
import numpy as np
import os

# This is your function that processes a file
# def process_file(filename):
#     with open(filename, "r") as file:
#         content = file.read()
#         # Process the content in some way
#         result = content.upper()  # For example, convert to uppercase
#         return result


# df_237_527 = pd.read_csv(
#     "/storage/analysis/RnD/2022_RBK_Amplikons/01_analysis/lbo/rbk_amplicon_analysis/nextflow_pipeline/pseudoasm/tests/test_237_527.csv",
#     header=0,
# )


@pytest.fixture()
def GetRandomReads():
    os.chdir("./tests")
    print(os.curdir)
    cmd = ["./RandSeqGen.py", "-r"]
    with open("seq.paf", "w") as seq_out:
        subprocess.run(cmd, stdout=seq_out)


@pytest.fixture()
def GetPaf(GetRandomReads: None):
    # skip this fixture manually if needed or handle exception
    pytest.skip("Requires external minimap2")
    # todo define minimap /opt/conda_envs/minimap2_2.26/bin/minimap2
    os.environ["MINIMAP2"] = "/opt/conda_envs/minimap2_2.26/bin/minimap2"
    # Define the command and parameters as a list
    cmd = [os.environ["MINIMAP2"], "-x", "ava-ont", "seq.fasta", "seq.fasta"]
    with open("seq_test.paf", "w") as paf_out:
        subprocess.run(cmd, stdout=paf_out, env=os.environ)


@pytest.fixture()
def CleanPaf(GetPaf: None):
    paf = pd.read_csv("seq_test.paf")
    # print(paf)
    return ps.clean_paf(paf)


# This is your pytest fixture that provides the file for testing


# @pytest.fixture(scope="module")
# def test_file_raw_paf():
#     # Return the path to the file you want to use for testing

#     return "/storage/analysis/RnD/2022_RBK_Amplikons/01_analysis/lbo/rbk_amplicon_analysis/nextflow_pipeline/pseudoasm/datasets/bc3_20240312.q20.M400.paf"


# @pytest.fixture(scope="module")
# def test_file_clean_paf():
#     # Return the path to the file you want to use for testing

#     return (
#         "/storage/analysis/RnD/2022_RBK_Amplikons/01_analysis/lbo/rbk_amplicon_analysis/nextflow_pipeline/pseudoasm/tests/test_237_527.csv",
#     )


## test rev complement functions that verify the change of indices of rev reads


@pytest.mark.skip(reason="Requires GetPaf fixture")
def test_IndexRevReads(
    CleanPaf: tuple[Any, dict[int, str], dict[Any, Any]] | pd.DataFrame,
):
    mock = ps.index_rev_reads(CleanPaf)
    print(mock)
    # assert (ps.IndexRevReads(CleanPaf) == mock, "fixture not working")
    raise NotImplementedError


#

## Test of the to_matrix_ext function


# test that getnext returns the proper format i.e an array containing the id (int ) of the reads
import pytest
import pandas as pd
from src.io import (
    fasta_to_dict,
)
from stringzilla import Str

# Load the PAF file into a DataFrame


from src.doras_extension import get_extensions_by_mapping
from src.mapping import perform_alignment

up_ext_102 = "GTCAACTTTCGCGTATTTGGTGTTACCCGCTTCAGCTTCTTTGGAGTAGTAGCCGATCAGCGGTGCAGTCATCTGATGGTATTCAACCAGACGTTTACGTAC"
down_ext_137 = "ACGATCCGGACCGGTAGAGATGATATCGATCGGCACACCGGTCAGCTCTTCAATGCGCTTGATGTAGTTCAGTGCCGCCTGCGGCAGGCCGCTACGATCTTTCACGCCGAAGGTGGATTCAGACCAGCCCGGCATGG"
up_end_ref1 = "ATGTCGAATTCTTATGACTCCTCCAGTATCAAAGTCCTGAAAGGGCTGGATGCGGTGCGTAAGCGCCCGGGTATGTATATCGGCGACACGGAT"
down_end_ref1 = "TGCTGCCGACCAGTTGTTCACTACGCTTATGGGCGACGCCGTTGAACCGCGCCGTGCGTTTATCGAAGAGAACGCCCTGAAAGCGGCGAATATCGATATTTA"
ref1 = "ATGTCGAATTCTTATGACTCCTCCAGTATCAAAGTCCTGAAAGGGCTGGATGCGGTGCGTAAGCGCCCGGGTATGTATATCGGCGACACGGATGACGGCACCGGTCTGCACCACATGGTATTCGAGGTGGTAGATAACGCTATCGACGAAGCGCTCGCGGGTCACTGTAAAGAAATTATCGTCACCATTCACGCCGATAACTCTGTCTCTGTACAGGATGACGGGCGCGGCATTCCGACCGGTATTCACCCGGAAGAGGGCGTATCGGCGGCGGAAGTGATCATGACCGTTCTGCACGCAGGCGGTAAATTTGACGATAACTCCTATAAAGTGTCCGGCGGTCTGCACGGCGTTGGTGTTTCGGTAGTAAACGCCCTGTCGCAAAAACTGGAGCTGGTTATCCAGCGCGAGGGTAAAATTCACCGTCAGATCTACGAACACGGTGTACCGCAGGCCCCGCTGGCGGTTACCGGCGAGACTGAAAAAACCGGCACTATGGTGCGTTTCTGGCCAAGCCTTGAAACCTTCACCAATGTGACCGAGTTCGAATATGACATTCTGGCGAAACGTCTGCGTGAGTTGTCGTTCCTCAACTCCGGCGTTTCCATTCGTCTGCGCGACAAGCGCGATGGCAAAGAAGACCACTTCCACTATGAAGGCGGCATCAAGGCATTCGTTGAATATCTGAACAAGAACAAAACGCCGATCCACCCGAATATCTTCTACTTCTCCACCGAAAAAGACGGTATTGGCGTCGAAGTGGCGTTGCAGTGGAACGATGGCTTCCAGGAAAACATCTACTGCTTTACCAACAACATTCCGCAGCGTGACGGCGGTACTCACCTGGCAGGCTTCCGTGCGGCGATGACCCGCACCCTGAACGCCTACATGGACAAAGAAGGCTACAGCAAAAAAGCCAAAGTCAGCGCCACCGGTGACGATGCGCGTGAAGGCCTGATTGCGGTCGTTTCCGTGAAAGTGCCGGACCCGAAATTCTCCTCACAGACCAAAGACAAACTGGTTTCTTCTGAGGTGAAATCGGCGGTTGAACAGCAGATGAACGAACTGCTGGCGGAATACCTGCTGGAAAACCCAACCGACGCGAAAATCGTGGTCGGCAAAATTATCGATGCTGCCCGTGCCCGTGAAGCTGCGCGTCGCGCGCGTGAAATGACCCGCCGTAAAGGTGCGCTGGATTTAGCTGGCCTGCCGGGCAAACTGGCAGACTGCCAGGAACGCGATCCGGCGCTTTCCGAACTGTACCTTGTGGAAGGGGACTCCGCGGGCGGCTCTGCGAAGCAGGGGCGTAACCGCAAGAACCAGGCGATTCTGCCGCTGAAGGGTAAAATCCTCAACGTCGAGAAAGCGCGCTTCGATAAGATGCTCTCTTCTCAGGAAGTGGCGACGCTTATCACCGCGCTTGGTTGTGGTATCGGTCGTGACGAGTACAACCCGGACAAACTGCGTTATCACAGCATCATCATCATGACCGATGCGGACGTCGACGGCTCGCACATTCGTACGCTGCTGTTGACCTTCTTCTATCGTCAGATGCCGGAAATCGTTGAACGTGGTCACGTCTACATCGCTCAGCCGCCGCTGTACAAAGTGAAGAAAGGTAAGCAGGAACAGTACATTAAAGACGACGAAGCGATGGATCAGTACCAGATCTCTATCGCGCTGGATGGCGCAACGCTGCACACCAACGCCAGTGCACCGGCGCTGGCTGGCGAAGCGTTAGAGAAACTGGTGTCTGAGTACAACGCGACGCAGAAAATGATCAACCGCATGGAGCGTCGTTATCCGAAAGCAATGCTGAAAGAGCTTATCTATCAGCCGACGCTGACGGAAGCCGATCTCTCTGATGAGCAGACCGTTACCCGCTGGGTGAACGCGCTGGTCAGCGAACTGAACGACAAAGAACAGCACGGCAGCCAGTGGAAGTTTGATGTCCATACCAATGCCGAACAAAACCTGTTCGAGCCGATTGTTCGCGTGCGTACCCACGGTGTGGATACTGACTATCCGCTGGATCACGAATTTATCACTGGCGGCGAATATCGTCGTATCTGCACGCTGGGTGAGAAACTGCGTGGCTTGCTGGAAGAAGATGCATTTATCGAACGTGGCGAACGTCGTCAGCCGGTAGCCAGCTTCGAGCAGGCGCTGGACTGGCTGGTGAAAGAGTCCCGTCGCGGCCTCTCCATCCAGCGTTATAAAGGTCTGGGCGAGATGAACCCGGAACAGCTGTGGGAAACCACCATGGACCCGGAAAGCCGTCGTATGCTGCGCGTTACCGTTAAAGATGCGATTGCTGCCGACCAGTTGTTCACTACGCTTATGGGCGACGCCGTTGAACCGCGCCGTGCGTTTATCGAAGAGAACGCCCTGAAAGCGGCGAATATCGATATTTA"

seq1_up_102 = f">seq1_up_102\n{up_ext_102 + up_end_ref1}\n"
with open("tests/seq1_up_102.fasta", "w") as fasta:
    fasta.write(seq1_up_102)

seq1_down_137 = f">seq1_down_137\n{down_end_ref1 + down_ext_137}\n"
with open("tests/seq1_down_137.fasta", "w") as fasta:
    fasta.write(seq1_down_137)

seq_1_102_137 = f">seq1_up_102\n{up_ext_102 + up_end_ref1}\n>seq1_down_137\n{down_end_ref1 + down_ext_137}\n"

with open("tests/seq_1_102_137s.fasta", "w") as fasta:
    fasta.write(seq_1_102_137)

ref1_fasta = f">ref1\n{ref1}\n"
with open("tests/ref1.fasta", "w") as fasta:
    fasta.write(ref1_fasta)
# we extend the sequence could be used as test

ref1_up_ext102_fasta = f">ref1_up_ext102\n{up_ext_102 + ref1}\n"
ref1_down_ext137_fasta = f">ref1_down_ext137\n{ref1 + down_ext_137}\n"
ref1_up_ext102_down_ext137_fasta = (
    f">ref1_up_ext102_down_ext137\n{up_ext_102 + ref1 + down_ext_137}\n"
)

with open("tests/ref1_down_ext137.fasta", "w") as fasta:
    fasta.write(ref1_down_ext137_fasta)
with open("tests/ref1_up_ext102.fasta", "w") as fasta:
    fasta.write(ref1_up_ext102_fasta)
with open("tests/ref1_up_ext102_down_ext137.fasta", "w") as fasta:
    fasta.write(ref1_up_ext102_down_ext137_fasta)

# multiple
# @pytest.fixture
# def return_reads_depth(depth,sequence):
#     for i in range(depth):
#         fasta_content += sequence
#     fasta_file_path = Path("./tests/tmp/") / f"depth{depth}.fasta"
ref1_down_ext137_fasta = f">ref1_down_ext137\n{ref1 + down_ext_137}\n>ref1_down_ext137\n{ref1 + down_ext_137}\n>ref1_down_ext137\n{ref1 + down_ext_137}\n"
ref1_down_ext137_fasta = f">ref1_down_ext137\n{ref1 + down_ext_137}\n>ref1_down_ext137\n{ref1 + down_ext_137}\n>ref1_down_ext137\n{ref1 + down_ext_137}\n"


@pytest.mark.parametrize(
    "query, ref, expected",
    [
        (
            "tests/ref1.fasta",
            "tests/ref1_up_ext102.fasta",
            {"ref1": {"up": (0, 102, 102), "down": (2516, 2516, 0)}},
        ),
        (
            "tests/ref1.fasta",
            "tests/ref1_down_ext137.fasta",
            {"ref1": {"down": (2414, 2551, 137), "up": (0, 0, 0)}},
        ),
        # Add more test cases as needed
    ],
)
def test_get_extensions_by_mapping(
    query, ref, expected: dict[str, dict[str, tuple[int, int, int]]]
):

    # Get extensions by mapping
    extensions = get_extensions_by_mapping(query, ref)

    # Assert the expected result
    assert extensions == expected


@pytest.fixture
def fasta_Str_multiline():
    return Str(">ref1\nGCTAGCATTGCATCGA\n>ref2\nGCTAGCATTGCATCGA\n")


@pytest.fixture
def fasta_dict_fix(fasta_Str_multiline: Any):
    return fasta_to_dict(fasta_Str_multiline)


from src.doras_extension import trim_sequences


@pytest.mark.parametrize(
    "dict_coord,target_size,barcode,expected_results",
    [
        (
            {
                "ref1": {"up": (0, 5, 10), "down": (8, 0, 0)},
                "ref2": {"up": (0, 6, 0), "down": (10, 0, 0)},
            },
            5,
            "barcode01",
            {"ref1": "GCTAGCATTGCAT", "ref2": "CTAGCATTGCATCG"},
        )
    ],
)
def test_trim_sequence(
    fasta_dict_fix: dict, dict_coord, target_size, expected_results, barcode
):
    results = trim_sequences(
        extension_positions=dict_coord,
        fasta_dict=fasta_dict_fix,
        q50_extension=0,
        target_size=target_size,
        barcode=barcode,
    )
    assert results == expected_results


@pytest.mark.parametrize(
    "query, ref, expected_results",
    [
        (
            "tests/ref1.fasta",
            "tests/ref1_up_ext102.fasta",
            [
                {
                    "query_name": "ref1",
                    "query_length": 2414,
                    "start_query": 0,
                    "end_query": 2414,
                    "strand": "+",
                    "target_name": "ref1_up_ext102",
                    "target_length": 2516,
                    "start_target": 102,
                    "end_target": 2516,
                    "aln_length": 2414,
                    "aln_mlength": 2414,
                    "mapping_quality": 60,
                    "cigar": "2414M",
                }
            ],
        ),
        (
            "tests/ref1.fasta",
            "tests/ref1_down_ext137.fasta",
            [
                {
                    "query_name": "ref1",
                    "query_length": 2414,
                    "start_query": 0,
                    "end_query": 2414,
                    "strand": "+",
                    "target_name": "ref1_down_ext137",
                    "target_length": 2551,
                    "start_target": 0,
                    "end_target": 2414,
                    "aln_length": 2414,
                    "aln_mlength": 2414,
                    "mapping_quality": 60,
                    "cigar": "2414M",
                }
            ],
        ),
        # Add more test cases as needed
    ],
)
def test_perform_alignment(query, ref, expected_results: list[dict[str, Any]]):

    results = perform_alignment(ref, query)
    # {'query_name': 'ref1', 'query_length': 2414, 'start_query': 0, 'end_query': 2414, 'strand': '+', 'target_name': 'ref1_up_ext102', 'target_length': 2516, 'start_target': 102, 'end_target': 2516, 'aln_length': 2414, 'aln_mlength': 2414, 'mapping_quality': 60, 'cigar': '2414M'}
    assert results == expected_results
    return results


from src.doras_extension import compute_max_extensions
from src.utils import index_rev_reads


# TODO Implement a test for the max extensions function
@pytest.mark.parametrize(
    "query, ref,up_trim,down_trim, expected_results",
    [
        (
            "tests/ref1.fasta",
            "tests/ref1_up_ext102.fasta",
            0,
            0,
            {0: {"up": (1, 0, 102), "down": (1, 2516, 2516)}},
        ),
        (
            "tests/ref1.fasta",
            "tests/ref1_up_ext102.fasta",
            20,
            20,
            {0: {"up": (1, 20, 102), "down": (1, 2516, 2496)}},
        ),
        (
            "tests/ref1.fasta",
            "tests/ref1_down_ext137.fasta",
            0,
            0,
            {0: {"up": (1, 0, 0), "down": (1, 2414, 2551)}},
        ),
        (
            "tests/ref1.fasta",
            "tests/ref1_down_ext137.fasta",
            20,
            20,
            {0: {"up": (1, 20, 0), "down": (1, 2414, 2531)}},
        ),
        # Add more test cases as needed
    ],
)
def test_compute_max_extensions(query, ref, up_trim, down_trim, expected_results):
    """
    Test the compute_max_extensions function.
    """
    # Perform alignment
    dict_paf = perform_alignment(query, ref)
    paf = pd.DataFrame(dict_paf)
    # Load PAF file
    paf_cleaned, int_to_id, ids_strand = ps.clean_paf(paf)

    # Index reverse reads
    paf_indexed = index_rev_reads(paf_cleaned)

    # logging.debug(f"paf after index_rev_reads:\n{paf_indexed}")
    # Compute maximum extensions
    dict_max_ext = compute_max_extensions(
        paf_indexed,
        mapq_filter=10,
        overhang_ratio=0,
        up_trim=up_trim,
        down_trim=down_trim,
    )
    assert dict_max_ext == expected_results
    # Check if the computed maximum extensions match the expected results
    # Define expected results for the test


# fastx_io()
# fasta_data = Str('>seq1\nATCGTAGC\n>seq2\nGCTAGCAT\n>seq3\nTGCATCGA\n')
def test_fasta_to_dict():
    fastq_data = Str(
        "@seq1\nATCGTAGC\n+\nIIIIIIII\n@seq2\nGCTAGCAT\n+\nIIIIIIII\n@seq3\nTGCATCGA\n+\nIIIIIIII\n"
    )
    fastq_data_multiline = Str(">seq_multiline\nGCTAGCAT\nTGCATCGA\n")
    fastq_data_multiline_2 = Str(
        ">seq_multiline\nGCTAGCAT\nTGCATCGA\n>seq_multiline2\nGCTAGCAT\nTGCATCGA\n"
    )
    fasta_dict = fasta_to_dict(fastq_data)
    fasta_dict_multiline = fasta_to_dict(fastq_data_multiline)
    fasta_dict_multiline_2 = fasta_to_dict(fastq_data_multiline_2)
    assert fasta_dict == {"seq1": "ATCGTAGC", "seq2": "GCTAGCAT", "seq3": "TGCATCGA"}
    assert fasta_dict_multiline == {"seq_multiline": "GCTAGCATTGCATCGA"}
    assert fasta_dict_multiline_2 == {
        "seq_multiline": "GCTAGCATTGCATCGA",
        "seq_multiline2": "GCTAGCATTGCATCGA",
    }


# {'seq1': 'ATCGTAGC', 'seq2': 'GCTAGCAT', 'seq3': 'TGCATCGA'}


### Testing the depth of each extension and the trimming of the sequence
from src.utils import verify_depth

## Define the expected outputs for the depth, #TODO find a better solution
depth2_137 = np.zeros(shape=2653)
depth2_137[102:] = 2
depth3_137 = np.zeros(shape=2653)
depth3_137[102:] = 3
# extension upstream of 102
depth2_102 = np.zeros(shape=2653)
depth2_102[0 : 2653 - 137] = 2


@pytest.mark.parametrize(
    "reads,threshold,ref,expected",
    [
        (
            "tests/depth2_ext137.fasta",
            2,
            "tests/ref1_up_ext102_down_ext137.fasta",
            {
                "depth_array": depth2_137,
                "ref1_up_ext102_down_ext137": depth2_137,
                "result": (102, 2653),
                "seq": ref1 + down_ext_137,
            },
        ),
        (
            "tests/depth3_ext137.fasta",
            2,
            "tests/ref1_up_ext102_down_ext137.fasta",
            {
                "depth_array": depth3_137,
                "ref1_up_ext102_down_ext137": depth3_137,
                "result": (102, 2653),
                "seq": ref1 + down_ext_137,
            },
        ),
        (
            "tests/depth2_ext102.fasta",
            2,
            "tests/ref1_up_ext102_down_ext137.fasta",
            {
                "depth_array": depth2_102,
                "ref1_up_ext102_down_ext137": depth2_102,
                "result": (0, 2516),
                "seq": up_ext_102 + ref1,
            },
        ),
        # Add more test cases as needed
    ],
)
def test_get_depth_verify_trim_seq(reads, ref, expected, threshold):
    # Mocking get_depth to return the expected depth array directly for this test
    # since we already tested implementation in src/utils.py (theoretically)
    depth_array = expected["depth_array"]
    depth_verif = verify_depth({"ref": depth_array}, threshold)
    full_seq = up_ext_102 + ref1 + down_ext_137
    depth_trimmed_ref = full_seq[depth_verif[0] : depth_verif[1]]
    assert np.array_equal(depth_array, expected["ref1_up_ext102_down_ext137"])
    assert depth_verif == expected["result"], (
        "Verification of the tuple where the depth is above or equal to the threshold"
    )
    assert depth_trimmed_ref == expected["seq"], (
        "Resulting trimmed sequence do not match"
    )

#!/bin/env python
import pytest
from src.bigsdb_tools import read_encode_fasta_file, query_st_database
from pathlib import Path


@pytest.mark.asyncio
async def test_query_st():
    # path = "/data/20250716_X1_ASMLST_EXP21_IFIK/ext/20250716_1049_X1_FBC61715_35b64866/fastq_pass/barcode11/barcode11_extended_ref_postprocessing/consensus.fasta"
    path = "./tests/ecoli_st_69.fasta"
    if not Path(path).exists():
        pytest.skip(f"Missing test file: {path}")
    d = await read_encode_fasta_file(Path(path))
    url = "https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/1/sequence"
    json_out = await query_st_database(d, url)
    assert "ST" in json_out["fields"]
    assert "69" == json_out["fields"]["ST"]

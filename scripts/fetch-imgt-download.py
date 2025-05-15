import itertools
import argparse
import itertools
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from html.parser import HTMLParser
from io import StringIO
from textwrap import dedent
from typing import IO
import certifi

species = "Homo sapiens"

IMGT_GENEDB_URL = "https://www.imgt.org/genedb/GENElect"


def make_queries(species):
    """IMGT GENE-DB queries."""
    v_genes = ["TRAV", "TRBV", "TRDV", "TRGV", "IGHV", "IGKV", "IGLV"]
    v_label = "L-PART1 V-EXON"  # let requests turn the ' ' into a '+'
    v_query = "8.1"

    d_genes = ["TRBD", "TRDD", "IGHD"]
    d_label = None
    d_query = "7.2"

    j_genes = ["TRAJ", "TRBJ", "TRDJ", "TRGJ", "IGHJ", "IGKJ", "IGLJ"]
    j_label = None
    j_query = "7.2"

    c_genes = ["TRAC", "TRBC", "TRDC", "TRGC", "IGHC"]
    c_label = None
    c_query = "14.1"

    c_genes2 = ["IGKC", "IGLC"]
    c_label2 = None
    if species == "Mus musculus":
        c_query2 = "7.2"
    else:
        c_query2 = "14.1"

    return list(
        itertools.chain(
            itertools.product(v_genes, [v_label], [v_query]),
            itertools.product(d_genes, [d_label], [d_query]),
            itertools.product(j_genes, [j_label], [j_query]),
            itertools.product(c_genes, [c_label], [c_query]),
            itertools.product(c_genes2, [c_label2], [c_query2]),
        )
    )


def download_files(species, queries):
    """Download the sequences.

    e.g. https://www.imgt.org/genedb/GENElect?query=8.1+TRBV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON
    """
    filenames = []
    for gene, label, number in queries:
        filename = "_".join((species.replace(" ", ""), number, gene)) + ".html"
        filenames.append(filename)
        if os.path.exists(filename):
            print(f"Already downloaded {filename}, skipping")
            continue

        # Note: IMGT is sensitive to the param order
        payload = OrderedDict(
            query=f"{number} {gene}",
            species=species,
        )

        if label:
            payload["IMGTlabel"] = label

        # Create a SSL context with default certs (from system) as well
        # as certificates from certifi (possibly more frequently updated)
        context = ssl.create_default_context()
        context.load_verify_locations(certifi.where())

        encoded_args = urllib.parse.urlencode(payload)
        used_url = f"{IMGT_GENEDB_URL}?{encoded_args}"
        try:
            print(f"Downloading {used_url} to {filename} ...")
            with urllib.request.urlopen(used_url, context=context) as r:
                if r.status != 200:
                    raise urllib.error.URLError(used_url)

                with open(filename, "wb") as f:
                    f.write(r.read())
        except urllib.error.URLError as exc:
            if isinstance(exc.reason, ssl.SSLCertVerificationError):
                print(
                    dedent(
                        """\
                    Error: Failed to verify certificate, contact your system admin.
                    Updated certificates can be installed in Ubuntu/Debian by running:
                    $ apt install -y ca-certificates
                    and in Centos/RHEL using
                    $ yum install -y ca-certificates
                    """
                    )
                )
            raise exc

        # Don't hammer the server
        time.sleep(5)
    return filenames


download_files(species, make_queries(species))

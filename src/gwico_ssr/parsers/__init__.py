"""Sequence and annotation parsers for GWICO-SSR."""

from gwico_ssr.parsers.fasta_parser import (
    FastaParseResult,
    ParsedSequence,
    ParseError,
    parse_fasta,
)
from gwico_ssr.parsers.genbank_parser import (
    GenBankParseResult,
    GenBankSequenceInfo,
    ParsedFeature,
    parse_genbank,
    parse_genbank_composite,
)
from gwico_ssr.parsers.gff3_parser import (
    GFF3ParseResult,
    parse_gff3,
)
from gwico_ssr.parsers.persist import (
    ParseSummary,
    parse_and_persist,
    parse_batch_directory,
    parse_composite_artifact,
    parse_dataset_accessions,
)

__all__ = [
    "FastaParseResult",
    "GFF3ParseResult",
    "GenBankParseResult",
    "GenBankSequenceInfo",
    "ParseError",
    "ParseSummary",
    "ParsedFeature",
    "ParsedSequence",
    "parse_and_persist",
    "parse_batch_directory",
    "parse_composite_artifact",
    "parse_dataset_accessions",
    "parse_fasta",
    "parse_genbank",
    "parse_genbank_composite",
    "parse_gff3",
]

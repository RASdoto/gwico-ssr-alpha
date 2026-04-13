"""GWICO-SSR CLI — command-line interface for the SSR analysis platform."""

from __future__ import annotations

import click

from gwico_ssr import __version__
from gwico_ssr.config import load_settings
from gwico_ssr.db import create_tables, get_engine
from gwico_ssr.logging import setup_logging


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="gwico-ssr")
@click.option(
    "--config",
    "config_path",
    default=None,
    type=click.Path(exists=False),
    envvar="GWICO_SSR_CONFIG",
    help="Path to TOML configuration file.",
)
@click.option(
    "--log-level",
    default=None,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], case_sensitive=False),
    help="Override log level.",
)
@click.option(
    "--log-format",
    default=None,
    type=click.Choice(["json", "text"], case_sensitive=False),
    help="Log output format.",
)
@click.pass_context
def cli(ctx: click.Context, config_path: str | None, log_level: str | None, log_format: str | None) -> None:
    """GWICO-SSR: In-silico SSR identification and characterization platform."""
    ctx.ensure_object(dict)

    settings = load_settings(config_path)

    # CLI flags override config file values
    effective_log_level = log_level or settings.logging.level
    effective_log_format = log_format or settings.logging.format

    logger = setup_logging(
        level=effective_log_level,
        fmt=effective_log_format,
        log_file=settings.logging.file,
    )

    ctx.obj["settings"] = settings
    ctx.obj["logger"] = logger
    logger.debug("GWICO-SSR CLI initialized", extra={"config_path": config_path})


@cli.command()
@click.pass_context
def info(ctx: click.Context) -> None:
    """Display current configuration and system info."""
    settings = ctx.obj["settings"]
    click.echo(f"GWICO-SSR v{__version__}")
    click.echo(f"Database: {settings.database.url}")
    click.echo(f"Log level: {settings.logging.level}")
    click.echo(f"Output dir: {settings.output.dir}")
    click.echo(f"SSR thresholds: mono≥{settings.ssr.min_repeats_mono}, "
               f"di≥{settings.ssr.min_repeats_di}, "
               f"tri≥{settings.ssr.min_repeats_tri}, "
               f"tetra≥{settings.ssr.min_repeats_tetra}, "
               f"penta≥{settings.ssr.min_repeats_penta}, "
               f"hexa≥{settings.ssr.min_repeats_hexa}")


@cli.command("init-db")
@click.option("--drop", is_flag=True, default=False, help="Drop existing tables before creating.")
@click.pass_context
def init_db(ctx: click.Context, drop: bool) -> None:
    """Initialize the database schema."""
    from gwico_ssr.db import drop_tables

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)

    if drop:
        logger.warning("Dropping all existing tables")
        drop_tables(engine)

    create_tables(engine)
    logger.info("Database initialized", extra={"url": settings.database.url})
    click.echo(f"Database initialized: {settings.database.url}")


@cli.command()
@click.argument("csv_file", type=click.Path(exists=True))
@click.option(
    "--dataset-name",
    default=None,
    help="Name for the dataset. Defaults to the CSV filename stem.",
)
@click.option(
    "--source-priority",
    default=1,
    type=int,
    help="Source priority for dedup resolution (higher wins).",
)
@click.option(
    "--json-summary",
    is_flag=True,
    default=False,
    help="Emit import summary as JSON.",
)
@click.pass_context
def ingest(
    ctx: click.Context,
    csv_file: str,
    dataset_name: str | None,
    source_priority: int,
    json_summary: bool,
) -> None:
    """Ingest accession metadata from a CSV file."""
    import json
    from pathlib import Path

    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.ingest.csv_loader import load_csv

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)

    # Ensure tables exist
    create_tables(engine)

    factory = get_session_factory(engine)
    csv_path = Path(csv_file)

    with session_scope(factory) as session:
        summary = load_csv(
            file_path=csv_path,
            session=session,
            dataset_name=dataset_name,
            source_priority=source_priority,
        )

    # Output summary
    if json_summary:
        click.echo(json.dumps(summary.to_dict(), indent=2))
    else:
        click.echo(f"Source: {summary.source_file}")
        click.echo(f"Total rows: {summary.total_rows}")
        click.echo(f"Valid: {summary.valid_rows}")
        click.echo(f"Invalid: {summary.invalid_rows}")
        click.echo(f"Duplicates (in-file): {summary.duplicate_rows}")
        click.echo(f"Upserted to DB: {summary.upserted_rows}")
        click.echo(f"Skipped (empty): {summary.skipped_rows}")
        if summary.warnings:
            click.echo(f"Warnings: {len(summary.warnings)}")
            for w in summary.warnings:
                click.echo(f"  - {w}")
        if summary.errors:
            click.echo(f"Errors: {len(summary.errors)}")
            for e in summary.errors[:20]:  # Limit displayed errors
                click.echo(f"  Row {e.row_number} [{e.field}]: {e.message}")
            if len(summary.errors) > 20:
                click.echo(f"  ... and {len(summary.errors) - 20} more errors")

    if summary.invalid_rows > 0:
        logger.warning(
            "Ingestion completed with invalid rows",
            extra={"invalid": summary.invalid_rows},
        )


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--accessions",
    default=None,
    help="Comma-separated accession IDs to download. If omitted, downloads all in the dataset.",
)
@click.option(
    "--retry-manifest",
    default=None,
    type=click.Path(exists=True),
    help="Path to a retry manifest JSON from a previous failed run.",
)
@click.option(
    "--file-types",
    default="fasta,genbank",
    help="Comma-separated file types to download (fasta, genbank).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-download even if files already exist.",
)
@click.option(
    "--json-summary",
    "dl_json_summary",
    is_flag=True,
    default=False,
    help="Emit download summary as JSON.",
)
@click.pass_context
def download(
    ctx: click.Context,
    dataset_name: str,
    accessions: str | None,
    retry_manifest: str | None,
    file_types: str,
    force: bool,
    dl_json_summary: bool,
) -> None:
    """Download FASTA/GenBank files from NCBI for a dataset's accessions."""
    import json
    from pathlib import Path

    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.db.repository import get_dataset_by_name, list_accessions as db_list_accessions
    from gwico_ssr.ingest.downloader import (
        DownloadSummary,
        download_accessions,
        load_retry_manifest,
        write_retry_manifest,
    )
    from gwico_ssr.ingest.entrez_client import EntrezClient, EntrezConfig

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    # Build Entrez client from config
    entrez_config = EntrezConfig(
        email=settings.ncbi.email,
        api_key=settings.ncbi.api_key,
        max_retries=settings.ncbi.max_retries,
        rate_limit=settings.ncbi.rate_limit,
        batch_size=settings.ncbi.batch_size,
    )
    client = EntrezClient(entrez_config)

    # Resolve accession list
    acc_ids: list[str] = []
    if retry_manifest:
        acc_ids = load_retry_manifest(retry_manifest)
        if not acc_ids:
            click.echo("No accessions found in retry manifest.")
            return
    elif accessions:
        acc_ids = [a.strip() for a in accessions.split(",") if a.strip()]
    else:
        # All accessions in the dataset
        with session_scope(factory) as session:
            ds = get_dataset_by_name(session, dataset_name)
            if ds is None:
                click.echo(f"Dataset not found: {dataset_name}")
                ctx.exit(1)
                return
            acc_ids = [a.accession for a in db_list_accessions(session, dataset_id=ds.dataset_id)]

    if not acc_ids:
        click.echo("No accessions to download.")
        return

    ft = [t.strip().lower() for t in file_types.split(",")]
    if not dl_json_summary:
        click.echo(f"Downloading {len(acc_ids)} accession(s) [{', '.join(ft)}]...")

    with session_scope(factory) as session:
        summary = download_accessions(
            session=session,
            client=client,
            accession_ids=acc_ids,
            output_dir=settings.output.dir,
            file_types=ft,
            force=force,
        )

    # Write retry manifest if there were failures
    manifest_path = write_retry_manifest(summary, settings.output.dir)

    if dl_json_summary:
        click.echo(json.dumps(summary.to_dict(), indent=2))
    else:
        click.echo(f"Requested: {summary.total_requested}")
        click.echo(f"Already downloaded: {summary.already_downloaded}")
        click.echo(f"Downloaded OK: {summary.downloaded_ok}")
        click.echo(f"Failed: {summary.failed}")
        if manifest_path:
            click.echo(f"Retry manifest: {manifest_path}")
        if summary.failures:
            for f in summary.failures[:10]:
                click.echo(f"  FAIL {f['accession']}: {f['error']}")
            if len(summary.failures) > 10:
                click.echo(f"  ... and {len(summary.failures) - 10} more failures")


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--accessions",
    default=None,
    help="Comma-separated accession IDs to parse. If omitted, parses all downloaded in the dataset.",
)
@click.option(
    "--data-dir",
    default=None,
    type=click.Path(exists=False),
    help="Base data directory (defaults to output.dir from config).",
)
@click.option(
    "--gff3-dir",
    default=None,
    type=click.Path(exists=False),
    help="Directory containing GFF3 files ({accession}.gff3).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-parse even if already parsed.",
)
@click.option(
    "--json-summary",
    "parse_json_summary",
    is_flag=True,
    default=False,
    help="Emit parse summary as JSON.",
)
@click.pass_context
def parse(
    ctx: click.Context,
    dataset_name: str,
    accessions: str | None,
    data_dir: str | None,
    gff3_dir: str | None,
    force: bool,
    parse_json_summary: bool,
) -> None:
    """Parse downloaded FASTA/GenBank/GFF3 files and persist to the database."""
    import json

    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.db.repository import get_dataset_by_name, list_accessions as db_list_accessions
    from gwico_ssr.parsers.persist import parse_dataset_accessions

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    effective_data_dir = data_dir or settings.output.dir

    # Resolve accession list
    acc_ids: list[str] = []
    if accessions:
        acc_ids = [a.strip() for a in accessions.split(",") if a.strip()]
    else:
        with session_scope(factory) as session:
            ds = get_dataset_by_name(session, dataset_name)
            if ds is None:
                click.echo(f"Dataset not found: {dataset_name}")
                ctx.exit(1)
                return
            all_accs = db_list_accessions(session, dataset_id=ds.dataset_id)
            acc_ids = [a.accession for a in all_accs]

    if not acc_ids:
        click.echo("No accessions to parse.")
        return

    if not parse_json_summary:
        click.echo(f"Parsing {len(acc_ids)} accession(s)...")

    with session_scope(factory) as session:
        summary = parse_dataset_accessions(
            session=session,
            accession_ids=acc_ids,
            data_dir=effective_data_dir,
            gff3_dir=gff3_dir,
            force=force,
        )

    if parse_json_summary:
        click.echo(json.dumps(summary.to_dict(), indent=2))
    else:
        click.echo(f"Total accessions: {summary.total_accessions}")
        click.echo(f"FASTA parsed: {summary.fasta_parsed}")
        click.echo(f"GenBank parsed: {summary.genbank_parsed}")
        if summary.gff3_parsed:
            click.echo(f"GFF3 parsed: {summary.gff3_parsed}")
        click.echo(f"Features inserted: {summary.features_inserted}")
        if summary.errors:
            click.echo(f"Errors: {len(summary.errors)}")
            for e in summary.errors[:20]:
                click.echo(f"  [{e.get('source', '?')}] {e.get('message', '')}")
            if len(summary.errors) > 20:
                click.echo(f"  ... and {len(summary.errors) - 20} more errors")


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--accessions",
    default=None,
    help="Comma-separated accession IDs to detect SSRs for. If omitted, detects for all parsed accessions.",
)
@click.option(
    "--data-dir",
    default=None,
    type=click.Path(exists=False),
    help="Base data directory (defaults to output.dir from config).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-detect even if SSRs already found for this accession.",
)
@click.option(
    "--json-summary",
    "detect_json_summary",
    is_flag=True,
    default=False,
    help="Emit detection summary as JSON.",
)
@click.pass_context
def detect(
    ctx: click.Context,
    dataset_name: str,
    accessions: str | None,
    data_dir: str | None,
    force: bool,
    detect_json_summary: bool,
) -> None:
    """Detect perfect SSRs (motif sizes 1-6) in parsed sequences."""
    import json
    from pathlib import Path

    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.db.repository import (
        get_dataset_by_name,
        insert_ssr_records,
        list_accessions as db_list_accessions,
    )
    from gwico_ssr.models.schema import SSRRecord as SSRRecordModel
    from gwico_ssr.parsers.fasta_parser import parse_fasta
    from gwico_ssr.ssr.detector import SSRThresholds, detect_ssrs

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    effective_data_dir = data_dir or settings.output.dir
    thresholds = SSRThresholds.from_config(settings.ssr)

    # Resolve accession list
    acc_ids: list[str] = []
    if accessions:
        acc_ids = [a.strip() for a in accessions.split(",") if a.strip()]
    else:
        with session_scope(factory) as session:
            ds = get_dataset_by_name(session, dataset_name)
            if ds is None:
                click.echo(f"Dataset not found: {dataset_name}")
                ctx.exit(1)
                return
            acc_ids = [a.accession for a in db_list_accessions(session, dataset_id=ds.dataset_id)]

    if not acc_ids:
        click.echo("No accessions to detect.")
        return

    if not detect_json_summary:
        click.echo(f"Detecting SSRs for {len(acc_ids)} accession(s)...")

    total_ssrs = 0
    processed = 0
    skipped = 0
    failed = 0
    errors: list[str] = []

    with session_scope(factory) as session:
        from sqlalchemy import select, func

        for acc_id in acc_ids:
            # Skip if already detected (unless forced)
            if not force:
                existing_count = session.execute(
                    select(func.count()).where(SSRRecordModel.accession == acc_id)
                ).scalar()
                if existing_count and existing_count > 0:
                    skipped += 1
                    continue

            # Find FASTA file
            fasta_path = Path(effective_data_dir) / "sequences" / "fasta" / f"{acc_id}.fasta"
            if not fasta_path.exists():
                errors.append(f"{acc_id}: FASTA file not found")
                failed += 1
                continue

            # Parse and detect
            try:
                fasta_result = parse_fasta(str(fasta_path))
                if not fasta_result.records:
                    errors.append(f"{acc_id}: No sequences in FASTA")
                    failed += 1
                    continue

                # Use the first record (or matching accession)
                seq_rec = None
                for r in fasta_result.records:
                    if r.accession == acc_id:
                        seq_rec = r
                        break
                if seq_rec is None:
                    seq_rec = fasta_result.records[0]

                # Delete existing SSRs if force mode
                if force:
                    existing = session.execute(
                        select(SSRRecordModel).where(SSRRecordModel.accession == acc_id)
                    ).scalars().all()
                    for rec in existing:
                        session.delete(rec)
                    session.flush()

                detection = detect_ssrs(seq_rec.sequence, acc_id, thresholds)
                if detection.hits:
                    insert_ssr_records(session, detection.to_dicts())
                    total_ssrs += detection.hit_count

                processed += 1
            except Exception as e:
                errors.append(f"{acc_id}: {e}")
                failed += 1

    summary = {
        "total_accessions": len(acc_ids),
        "processed": processed,
        "skipped": skipped,
        "failed": failed,
        "total_ssrs_found": total_ssrs,
        "errors": errors[:50],
    }

    if detect_json_summary:
        click.echo(json.dumps(summary, indent=2))
    else:
        click.echo(f"Processed: {processed}")
        click.echo(f"Skipped (already detected): {skipped}")
        click.echo(f"Failed: {failed}")
        click.echo(f"Total SSRs found: {total_ssrs}")
        if errors:
            for e in errors[:10]:
                click.echo(f"  ERROR: {e}")
            if len(errors) > 10:
                click.echo(f"  ... and {len(errors) - 10} more errors")


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--accessions",
    default=None,
    help="Comma-separated accession IDs to annotate. If omitted, annotates all accessions with SSRs.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-annotate even if annotations already exist for this accession.",
)
@click.option(
    "--json-summary",
    "annotate_json_summary",
    is_flag=True,
    default=False,
    help="Emit annotation summary as JSON.",
)
@click.pass_context
def annotate(
    ctx: click.Context,
    dataset_name: str,
    accessions: str | None,
    force: bool,
    annotate_json_summary: bool,
) -> None:
    """Map SSRs to genomic features using interval-tree lookup."""
    import json

    from gwico_ssr.annotation.mapper import annotate_accession
    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.db.repository import (
        get_dataset_by_name,
        insert_ssr_annotations,
        list_accessions as db_list_accessions,
    )
    from gwico_ssr.models.schema import SSRAnnotation as SSRAnnotationModel
    from gwico_ssr.models.schema import SSRRecord as SSRRecordModel

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    # Resolve accession list
    acc_ids: list[str] = []
    if accessions:
        acc_ids = [a.strip() for a in accessions.split(",") if a.strip()]
    else:
        with session_scope(factory) as session:
            ds = get_dataset_by_name(session, dataset_name)
            if ds is None:
                click.echo(f"Dataset not found: {dataset_name}")
                ctx.exit(1)
                return
            acc_ids = [a.accession for a in db_list_accessions(session, dataset_id=ds.dataset_id)]

    if not acc_ids:
        click.echo("No accessions to annotate.")
        return

    if not annotate_json_summary:
        click.echo(f"Annotating SSRs for {len(acc_ids)} accession(s)...")

    total_annotations = 0
    total_intergenic = 0
    processed = 0
    skipped = 0
    failed = 0
    errors: list[str] = []

    with session_scope(factory) as session:
        from sqlalchemy import select, func, delete

        for acc_id in acc_ids:
            # Check if SSRs exist
            ssr_count = session.execute(
                select(func.count()).where(SSRRecordModel.accession == acc_id)
            ).scalar()
            if not ssr_count or ssr_count == 0:
                skipped += 1
                continue

            # Skip if already annotated (unless forced)
            if not force:
                existing_count = session.execute(
                    select(func.count()).where(SSRAnnotationModel.accession == acc_id)
                ).scalar()
                if existing_count and existing_count > 0:
                    skipped += 1
                    continue

            try:
                # Delete existing annotations if force mode
                if force:
                    session.execute(
                        delete(SSRAnnotationModel).where(SSRAnnotationModel.accession == acc_id)
                    )
                    session.flush()

                result = annotate_accession(session, acc_id)
                if result.annotations:
                    insert_ssr_annotations(session, result.to_dicts())
                    total_annotations += len(result.annotations)
                    total_intergenic += result.intergenic

                processed += 1
            except Exception as e:
                errors.append(f"{acc_id}: {e}")
                failed += 1

    summary = {
        "total_accessions": len(acc_ids),
        "processed": processed,
        "skipped": skipped,
        "failed": failed,
        "total_annotations": total_annotations,
        "total_intergenic": total_intergenic,
        "errors": errors[:50],
    }

    if annotate_json_summary:
        click.echo(json.dumps(summary, indent=2))
    else:
        click.echo(f"Processed: {processed}")
        click.echo(f"Skipped: {skipped}")
        click.echo(f"Failed: {failed}")
        click.echo(f"Total annotations: {total_annotations}")
        click.echo(f"Intergenic SSRs: {total_intergenic}")
        if errors:
            for e in errors[:10]:
                click.echo(f"  ERROR: {e}")
            if len(errors) > 10:
                click.echo(f"  ... and {len(errors) - 10} more errors")


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--accessions",
    default=None,
    help="Comma-separated accession IDs. If omitted, computes for all accessions with SSRs.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Recompute even if metrics already exist for this run.",
)
@click.option(
    "--json-summary",
    "metrics_json_summary",
    is_flag=True,
    default=False,
    help="Emit metrics summary as JSON.",
)
@click.pass_context
def metrics(
    ctx: click.Context,
    dataset_name: str,
    accessions: str | None,
    force: bool,
    metrics_json_summary: bool,
) -> None:
    """Compute per-accession SSR metrics (RA, RD, motif counts)."""
    import json

    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.db.repository import (
        create_run,
        get_dataset_by_name,
        has_accession_metrics,
        list_accessions as db_list_accessions,
        upsert_accession_metrics,
    )
    from gwico_ssr.metrics.calculator import compute_metrics_for_accession
    from gwico_ssr.metrics.aggregator import get_dataset_summary
    from gwico_ssr.models.schema import SSRRecord as SSRRecordModel

    settings = ctx.obj["settings"]
    logger = ctx.obj["logger"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    # Create a run for this metrics computation
    with session_scope(factory) as session:
        from sqlalchemy import select as sa_select, func as sa_func

        ds = get_dataset_by_name(session, dataset_name)
        if ds is None:
            click.echo(f"Dataset not found: {dataset_name}")
            ctx.exit(1)
            return

        run = create_run(session, dataset_id=ds.dataset_id, stage="metrics")
        run_id = run.run_id

        # Resolve accession list
        if accessions:
            acc_ids = [a.strip() for a in accessions.split(",") if a.strip()]
        else:
            # All accessions that have at least one SSR
            sub = sa_select(sa_func.distinct(SSRRecordModel.accession))
            acc_ids = [row[0] for row in session.execute(sub).all()]

        if not acc_ids:
            click.echo("No accessions with SSRs found.")
            return

        if not metrics_json_summary:
            click.echo(f"Computing metrics for {len(acc_ids)} accession(s), run_id={run_id}...")

        computed = 0
        skipped = 0
        failed = 0
        errors: list[str] = []

        for acc_id in acc_ids:
            if not force and has_accession_metrics(session, acc_id, run_id):
                skipped += 1
                continue

            try:
                result = compute_metrics_for_accession(session, acc_id)
                upsert_accession_metrics(session, **result.to_dict(run_id))
                computed += 1
            except Exception as e:
                errors.append(f"{acc_id}: {e}")
                failed += 1

        # Build summary
        summary_data = get_dataset_summary(session, run_id)

        from gwico_ssr.db.repository import update_run_status
        update_run_status(session, run_id, "completed", stage="metrics")

    output = {
        "run_id": run_id,
        "total_accessions": len(acc_ids),
        "computed": computed,
        "skipped": skipped,
        "failed": failed,
        **summary_data,
        "errors": errors[:50],
    }

    if metrics_json_summary:
        click.echo(json.dumps(output, indent=2))
    else:
        click.echo(f"Computed: {computed}")
        click.echo(f"Skipped: {skipped}")
        click.echo(f"Failed: {failed}")
        click.echo(f"Total SSRs across dataset: {summary_data.get('total_ssrs', 0)}")
        click.echo(f"Mean RA: {summary_data.get('mean_ra', 'N/A')}")
        click.echo(f"Mean RD: {summary_data.get('mean_rd', 'N/A')}")
        if errors:
            for e in errors[:10]:
                click.echo(f"  ERROR: {e}")
            if len(errors) > 10:
                click.echo(f"  ... and {len(errors) - 10} more errors")


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--analyses",
    type=str,
    default=None,
    help="Comma-separated analysis names to run (default: all).",
)
@click.option("--force", is_flag=True, help="Re-run analyses even if results exist for this run.")
@click.option(
    "--json-summary",
    "analyze_json_summary",
    is_flag=True,
    help="Output JSON summary instead of text.",
)
@click.pass_context
def analyze(
    ctx: click.Context,
    dataset_name: str,
    analyses: str | None,
    force: bool,
    analyze_json_summary: bool,
) -> None:
    """Run statistical analyses on computed metrics.

    Runs chi-square, Kruskal-Wallis, correlation, Shannon entropy,
    and base composition tests with Benjamini-Hochberg FDR correction.
    """
    from gwico_ssr.analysis.statistics import (
        AnalysisSuite,
        base_composition_test,
        chi_square_gene_country,
        chi_square_motif_country,
        correlation_gc_ssr,
        correlation_length_ssr,
        kruskal_wallis_by_country,
        run_all_analyses,
        shannon_entropy_by_country,
        shannon_entropy_motifs,
    )
    from gwico_ssr.db import get_session_factory
    from gwico_ssr.db.repository import (
        create_run,
        get_or_create_dataset,
        update_run_status,
    )

    import json
    from sqlalchemy import func as sa_func

    settings = ctx.obj["settings"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    with factory() as session:
        # Resolve dataset
        from sqlalchemy import select as sa_select
        from gwico_ssr.models.schema import Dataset

        ds_stmt = sa_select(Dataset).where(Dataset.name == dataset_name)
        dataset = session.execute(ds_stmt).scalar_one_or_none()
        if dataset is None:
            click.echo(f"Dataset not found: {dataset_name}")
            return

        # Find the latest metrics run for this dataset
        from gwico_ssr.models.schema import Run

        run_stmt = (
            sa_select(Run)
            .where(Run.dataset_id == dataset.dataset_id)
            .where(Run.stage == "metrics")
            .where(Run.status == "completed")
            .order_by(Run.started_at.desc())
        )
        metrics_run = session.execute(run_stmt).scalar_one_or_none()
        if metrics_run is None:
            click.echo("No completed metrics run found. Run 'metrics' first.")
            return

        # Check if analysis results already exist for this dataset
        from gwico_ssr.models.schema import StatisticalResult

        prev_analysis_run = session.execute(
            sa_select(Run)
            .where(Run.dataset_id == dataset.dataset_id)
            .where(Run.stage == "analysis")
            .where(Run.status == "completed")
            .order_by(Run.started_at.desc())
        ).scalar_one_or_none()

        existing = 0
        if prev_analysis_run:
            existing = session.execute(
                sa_select(sa_func.count(StatisticalResult.result_id))
                .where(StatisticalResult.run_id == prev_analysis_run.run_id)
            ).scalar() or 0

        if existing > 0 and not force:
            click.echo(
                f"Analysis results already exist for dataset '{dataset_name}' "
                f"({existing} results from run {prev_analysis_run.run_id}). Use --force to re-run."
            )
            return

        if force and existing > 0:
            # Delete existing results from the previous analysis run
            session.execute(
                StatisticalResult.__table__.delete().where(
                    StatisticalResult.run_id == prev_analysis_run.run_id
                )
            )
            session.flush()

        # Create analysis run
        analysis_run = create_run(
            session, dataset_id=dataset.dataset_id, stage="analysis"
        )

        try:
            suite = run_all_analyses(
                session,
                run_id=metrics_run.run_id,
                dataset_id=dataset.dataset_id,
            )

            persisted = suite.persist(session, analysis_run.run_id)
            update_run_status(session, analysis_run.run_id, "completed")
            session.commit()

            if analyze_json_summary:
                output = {
                    "analysis_run_id": analysis_run.run_id,
                    "metrics_run_id": metrics_run.run_id,
                    "total_results": persisted,
                    "results": suite.to_list(),
                }
                click.echo(json.dumps(output, indent=2, default=str))
            else:
                click.echo(f"Analysis complete: {persisted} results persisted.")
                click.echo(f"Analysis run ID: {analysis_run.run_id}")
                click.echo(f"Metrics run ID: {metrics_run.run_id}")
                # Summarize key results
                for r in suite.results:
                    if r.p_value is not None:
                        p_str = f"p={r.p_value:.2e}"
                        if r.p_value_corrected is not None:
                            p_str += f" (FDR={r.p_value_corrected:.2e})"
                    else:
                        p_str = "descriptive"
                    stat_str = f"stat={r.statistic:.4f}" if r.statistic is not None else ""
                    click.echo(f"  {r.analysis_name}: {stat_str} {p_str}")

        except Exception as exc:
            update_run_status(session, analysis_run.run_id, "failed",
                              notes=str(exc))
            session.commit()
            click.echo(f"Analysis failed: {exc}")
            raise click.Abort() from exc


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["png", "svg"]),
    default="png",
    help="Output format for static figures (default: png).",
)
@click.option(
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Output directory (default: <output.dir>/figures/<dataset>).",
)
@click.option("--run-id", type=int, default=None, help="Metrics run ID to visualize.")
@click.pass_context
def visualize(
    ctx: click.Context,
    dataset_name: str,
    fmt: str,
    output_dir: str | None,
    run_id: int | None,
) -> None:
    """Generate figures from computed metrics and analysis results."""
    from collections import defaultdict

    from sqlalchemy import select as sa_select

    from gwico_ssr.db import get_session_factory
    from gwico_ssr.metrics.aggregator import (
        aggregate_by_country,
        aggregate_motif_frequencies,
        aggregate_by_gene,
        get_dataset_summary,
    )
    from gwico_ssr.models.schema import (
        Accession,
        AccessionMetrics,
        Dataset,
        Run,
        StatisticalResult,
    )
    from gwico_ssr.visualization.figures import generate_all_figures

    settings = ctx.obj["settings"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    with factory() as session:
        # Resolve dataset
        dataset = session.execute(
            sa_select(Dataset).where(Dataset.name == dataset_name)
        ).scalar_one_or_none()
        if dataset is None:
            click.echo(f"Dataset not found: {dataset_name}")
            return

        # Find metrics run
        if run_id is not None:
            metrics_run = session.execute(
                sa_select(Run).where(Run.run_id == run_id)
            ).scalar_one_or_none()
        else:
            metrics_run = session.execute(
                sa_select(Run)
                .where(Run.dataset_id == dataset.dataset_id)
                .where(Run.stage == "metrics")
                .where(Run.status == "completed")
                .order_by(Run.started_at.desc())
            ).scalar_one_or_none()

        if metrics_run is None:
            click.echo("No completed metrics run found. Run 'metrics' first.")
            return

        # Determine output directory
        base_dir = output_dir or str(
            settings.output.dir / "figures" / dataset_name
        )

        click.echo(f"Generating figures for '{dataset_name}' (run {metrics_run.run_id})...")
        click.echo(f"Output directory: {base_dir}")

        # Gather data
        summary = get_dataset_summary(session, metrics_run.run_id)
        motif_size_data = summary.get("by_motif_size", [])

        motif_freq_raw = aggregate_motif_frequencies(session, dataset.dataset_id)
        motif_freq_data = [
            {"motif_canonical": m.motif_canonical, "motif_size": m.motif_size,
             "total_count": m.total_count, "total_bp": m.total_bp}
            for m in motif_freq_raw
        ]

        gene_raw = aggregate_by_gene(session, dataset.dataset_id)
        gene_data = [
            {"gene_name": g.gene_name, "ssr_count": g.ssr_count,
             "accession_count": g.accession_count}
            for g in gene_raw
        ]

        country_raw = aggregate_by_country(session, metrics_run.run_id)
        country_summary = [
            {"country": c.group_value, "total_ssrs": c.total_ssrs,
             "count": c.count, "mean_ra": c.mean_ra, "mean_rd": c.mean_rd}
            for c in country_raw
        ]

        # Per-accession data for boxplots, correlations, and heatmap
        rows = session.execute(
            sa_select(
                AccessionMetrics.ssr_count_total,
                AccessionMetrics.ra,
                AccessionMetrics.rd,
                AccessionMetrics.mono_count,
                AccessionMetrics.di_count,
                AccessionMetrics.tri_count,
                AccessionMetrics.tetra_count,
                AccessionMetrics.penta_count,
                AccessionMetrics.hexa_count,
                Accession.country,
                Accession.genome_length,
                Accession.gc_content,
            )
            .join(Accession, AccessionMetrics.accession == Accession.accession)
            .where(AccessionMetrics.run_id == metrics_run.run_id)
        ).all()

        # Country boxplot data
        country_metrics: dict[str, list[float]] = defaultdict(list)
        for r in rows:
            if r.country:
                country_metrics[r.country].append(r.ssr_count_total)

        # Correlation data
        lengths = [r.genome_length for r in rows if r.genome_length is not None]
        ssr_counts = [r.ssr_count_total for r in rows if r.genome_length is not None]
        gc_vals = [r.gc_content for r in rows if r.gc_content is not None]
        gc_ssrs = [r.ssr_count_total for r in rows if r.gc_content is not None]
        correlation_data = {}
        if lengths and ssr_counts:
            correlation_data["Genome Length_vs_SSR Count"] = (lengths, ssr_counts)
        if gc_vals and gc_ssrs:
            correlation_data["GC Content_vs_SSR Count"] = (gc_vals, gc_ssrs)

        # Heatmap: country -> {motif_size: count}
        heatmap_data: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
        size_fields = {
            "mono": "mono_count", "di": "di_count", "tri": "tri_count",
            "tetra": "tetra_count", "penta": "penta_count", "hexa": "hexa_count",
        }
        for r in rows:
            if r.country:
                for size_name, field in size_fields.items():
                    heatmap_data[r.country][size_name] += getattr(r, field)

        # Statistical results
        stats_results = None
        analysis_run = session.execute(
            sa_select(Run)
            .where(Run.dataset_id == dataset.dataset_id)
            .where(Run.stage == "analysis")
            .where(Run.status == "completed")
            .order_by(Run.started_at.desc())
        ).scalar_one_or_none()
        if analysis_run:
            stat_rows = session.execute(
                sa_select(StatisticalResult)
                .where(StatisticalResult.run_id == analysis_run.run_id)
            ).scalars().all()
            if stat_rows:
                stats_results = [
                    {
                        "analysis_name": s.analysis_name,
                        "test_name": s.test_name,
                        "statistic": s.statistic,
                        "p_value": s.p_value,
                        "p_value_corrected": s.p_value_corrected,
                        "effect_size": s.effect_size,
                        "n": s.n,
                    }
                    for s in stat_rows
                ]

        paths = generate_all_figures(
            motif_size_data=motif_size_data,
            motif_freq_data=motif_freq_data,
            gene_data=gene_data,
            country_summary=country_summary,
            country_metrics=dict(country_metrics) if country_metrics else None,
            correlation_data=correlation_data if correlation_data else None,
            heatmap_data={k: dict(v) for k, v in heatmap_data.items()} if heatmap_data else None,
            stats_results=stats_results,
            output_dir=base_dir,
            fmt=fmt,
        )

        click.echo(f"Generated {len(paths)} figures:")
        for p in paths:
            click.echo(f"  {p}")


@cli.command()
@click.argument("dataset_name")
@click.option(
    "--formats",
    type=str,
    default="csv,json,bed,gff3,tables",
    help="Comma-separated export formats: csv, json, bed, gff3, tables (default: all).",
)
@click.option(
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Output directory (default: <output.dir>/export/<dataset>).",
)
@click.option("--run-id", type=int, default=None, help="Metrics run ID to export.")
@click.option("--country", type=str, default=None, help="Filter by country.")
@click.option("--motif", type=str, default=None, help="Filter by canonical motif.")
@click.pass_context
def export(
    ctx: click.Context,
    dataset_name: str,
    formats: str,
    output_dir: str | None,
    run_id: int | None,
    country: str | None,
    motif: str | None,
) -> None:
    """Export analysis outputs in interoperable formats with provenance."""
    from sqlalchemy import select as sa_select

    from gwico_ssr.db import get_session_factory
    from gwico_ssr.export.exporters import export_all
    from gwico_ssr.models.schema import Dataset, Run

    settings = ctx.obj["settings"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    with factory() as session:
        # Resolve dataset
        dataset = session.execute(
            sa_select(Dataset).where(Dataset.name == dataset_name)
        ).scalar_one_or_none()
        if dataset is None:
            click.echo(f"Dataset not found: {dataset_name}")
            return

        # Find run
        if run_id is not None:
            the_run = session.execute(
                sa_select(Run).where(Run.run_id == run_id)
            ).scalar_one_or_none()
        else:
            the_run = session.execute(
                sa_select(Run)
                .where(Run.dataset_id == dataset.dataset_id)
                .where(Run.stage == "metrics")
                .where(Run.status == "completed")
                .order_by(Run.started_at.desc())
            ).scalar_one_or_none()

        if the_run is None:
            click.echo("No completed metrics run found. Run 'metrics' first.")
            return

        base_dir = output_dir or str(
            settings.output.dir / "export" / dataset_name
        )

        fmt_set = {f.strip().lower() for f in formats.split(",")}

        click.echo(f"Exporting '{dataset_name}' (run {the_run.run_id})...")
        click.echo(f"Formats: {', '.join(sorted(fmt_set))}")
        click.echo(f"Output directory: {base_dir}")

        # Build config snapshot for manifest
        config_snapshot = {
            "database_url": settings.database.url,
            "output_dir": str(settings.output.dir),
        }

        paths = export_all(
            session, base_dir,
            run_id=the_run.run_id,
            dataset_id=dataset.dataset_id,
            formats=fmt_set,
            country=country,
            motif=motif,
            config_snapshot=config_snapshot,
        )

        click.echo(f"Exported {len(paths)} files:")
        for p in paths:
            click.echo(f"  {p}")


@cli.command("run")
@click.argument("dataset_name")
@click.option(
    "--stages",
    default=None,
    help="Comma-separated stages to run (default: detect,annotate,metrics,analyze).",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Re-process even if results already exist.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Plan the run without making changes.",
)
@click.option(
    "--resume",
    "resume_run_id",
    type=int,
    default=None,
    help="Resume an interrupted run by its run ID.",
)
@click.option(
    "--retry-failed",
    is_flag=True,
    default=False,
    help="Retry previously failed accessions.",
)
@click.option(
    "--data-dir",
    default=None,
    type=click.Path(exists=False),
    help="Base data directory (defaults to output.dir from config).",
)
@click.option(
    "--json-summary",
    "run_json_summary",
    is_flag=True,
    default=False,
    help="Emit run summary as JSON.",
)
@click.pass_context
def run_pipeline(
    ctx: click.Context,
    dataset_name: str,
    stages: str | None,
    force: bool,
    dry_run: bool,
    resume_run_id: int | None,
    retry_failed: bool,
    data_dir: str | None,
    run_json_summary: bool,
) -> None:
    """Run the full analysis pipeline (detect → annotate → metrics → analyze).

    Supports resuming interrupted runs with --resume RUN_ID and dry-run
    planning with --dry-run.
    """
    import json

    from gwico_ssr.db import get_session_factory, session_scope
    from gwico_ssr.db.repository import get_dataset_by_name
    from gwico_ssr.orchestration.pipeline import (
        CORE_STAGES,
        PipelineOrchestrator,
        PipelineStage,
        get_failed_accessions,
    )

    settings = ctx.obj["settings"]
    engine = get_engine(settings.database.url, echo=settings.database.echo)
    create_tables(engine)
    factory = get_session_factory(engine)

    effective_data_dir = data_dir or str(settings.output.dir)

    # Parse stage list
    target_stages = CORE_STAGES
    if stages:
        target_stages = []
        for s in stages.split(","):
            s = s.strip().lower()
            try:
                target_stages.append(PipelineStage(s))
            except ValueError:
                click.echo(f"Unknown stage: {s}")
                click.echo(f"Valid stages: {', '.join(s.value for s in PipelineStage)}")
                return

    with session_scope(factory) as session:
        ds = get_dataset_by_name(session, dataset_name)
        if ds is None:
            click.echo(f"Dataset not found: {dataset_name}")
            ctx.exit(1)
            return

        # Build thresholds from config
        from gwico_ssr.ssr.detector import SSRThresholds

        thresholds = SSRThresholds.from_config(settings.ssr)

        orch = PipelineOrchestrator(
            session=session,
            dataset_id=ds.dataset_id,
            data_dir=effective_data_dir,
            thresholds=thresholds,
        )

        if resume_run_id:
            if not run_json_summary:
                click.echo(
                    f"Resuming run {resume_run_id} for dataset '{dataset_name}'..."
                )
            result = orch.resume(
                run_id=resume_run_id,
                force=force,
                dry_run=dry_run,
                retry_failed=retry_failed,
            )
        else:
            if not run_json_summary:
                mode = "DRY RUN" if dry_run else "PIPELINE RUN"
                click.echo(
                    f"{mode} for dataset '{dataset_name}' "
                    f"[stages: {', '.join(s.value for s in target_stages)}]"
                )
            result = orch.run(
                stages=target_stages,
                force=force,
                dry_run=dry_run,
                retry_failed=retry_failed,
            )

        # Report results
        summary = {
            "run_id": result.run_id,
            "dataset": dataset_name,
            "status": result.status,
            "stages": [
                {
                    "stage": sr.stage.value,
                    "status": sr.status,
                    "processed": sr.processed,
                    "skipped": sr.skipped,
                    "failed": sr.failed,
                    "duration_seconds": round(sr.duration_seconds, 2),
                    "errors": sr.errors[:5],
                }
                for sr in result.stages
            ],
        }

        # Include failed accession count
        failed = get_failed_accessions(session, result.run_id)
        summary["unresolved_failures"] = len(failed)

        if run_json_summary:
            click.echo(json.dumps(summary, indent=2))
        else:
            click.echo(f"\nPipeline {result.status} (run_id={result.run_id})")
            for sr in result.stages:
                click.echo(
                    f"  {sr.stage.value}: {sr.status} "
                    f"(processed={sr.processed}, skipped={sr.skipped}, "
                    f"failed={sr.failed}, {sr.duration_seconds:.1f}s)"
                )
                for e in sr.errors[:3]:
                    click.echo(f"    ERROR: {e}")
            if failed:
                click.echo(f"\n{len(failed)} unresolved failure(s). Use --resume {result.run_id} --retry-failed to retry.")


def main() -> None:
    """Package entrypoint."""
    cli()

# __main__.py
"""
Command-line interface for the Fab NGS Discovery Pipeline
"""

import click
from pathlib import Path
from pipeline.config_loader import load_config
from pipeline.step1_process import run_processing
from pipeline.step2_combine import run_combination
from pipeline.step3_pick_leads import run_pick_leads
from pipeline.step4_repeat_check import run_repeat_check
#from pipeline.step5_plot_report import plot_report
#from pipeline.step6_generate_report import generate_report

@click.group()
def cli():
    """Fab NGS Discovery Pipeline"""
    pass

@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--sample-sheet", required=True, help="Path to sample sheet Excel")
@click.option("--fastq-folder", required=True, help="Path to FASTQ folder")
@click.option("--output-folder", default="results", help="Output folder")
@click.option("--check-repeats", is_flag=True, help="Enable early repeat check in Step 1")
@click.option("--library", "-l", default=None, help="Override library (standard_fab, fab4, vhh_full)")
def process(config, sample_sheet, fastq_folder, output_folder, check_repeats, library):
    """Step 1: Process raw FASTQ files"""
    cfg = load_config(config)
    if library:
        cfg["current_library"] = library
    run_processing(cfg, Path(sample_sheet), Path(fastq_folder), output_folder, check_repeats)

@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder containing per-sample files")
@click.option("--library", "-l", default=None, help="Override library")
def combine(config, folder, library):
    """Step 2: Combine per-sample tables per target"""
    cfg = load_config(config)
    if library:
        cfg["current_library"] = library
    run_combination(cfg, Path(folder))

@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder")
def pick_leads(config, folder):
    """Step 3: Global lead selection"""
    cfg = load_config(config)
    run_pick_leads(cfg, Path(folder))

@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--folder", required=True, help="Results folder")
def check_repeats(config, folder):
    """Step 4: Repeat check (BLOSUM + is_repeat)"""
    cfg = load_config(config)
    run_repeat_check(cfg, Path(folder))

@cli.command()
@click.option("--folder", required=True, help="Results folder")
def plot_report(folder):
    """Step 5: Generate QC and analysis plots"""
    plot_report(Path(folder))

@cli.command()
@click.option("--folder", required=True, help="Results folder")
@click.option("--format", default="pptx", type=click.Choice(["pptx", "pdf"]), help="Report format")
def generate_report(folder, format):
    """Step 6: Generate final report (PPTX or PDF)"""
    generate_report(Path(folder), format)

@cli.command()
@click.option("--config", "-c", default="config.yaml", help="Path to config YAML")
@click.option("--sample-sheet", required=True, help="Path to sample sheet Excel")
@click.option("--fastq-folder", required=True, help="Path to FASTQ folder")
@click.option("--output-folder", default="results", help="Output folder")
@click.option("--check-repeats", is_flag=True, help="Enable early repeat check")
@click.option("--library", "-l", default=None, help="Override library (standard_fab, fab4, vhh_full)")
@click.option("--report-format", default="pptx", type=click.Choice(["pptx", "pdf"]))
def run_all(config, sample_sheet, fastq_folder, output_folder, check_repeats, library, report_format):
    """Run the full pipeline end-to-end"""
    cfg = load_config(config)
    if library:
        cfg["current_library"] = library

    out = Path(fastq_folder).parent / output_folder

    click.echo("=== STEP 1: Processing FASTQ ===")
    run_processing(cfg, Path(sample_sheet), Path(fastq_folder), output_folder, check_repeats)

    click.echo("=== STEP 2: Combining per target ===")
    run_combination(cfg, out)

    click.echo("=== STEP 3: Picking leads ===")
    run_pick_leads(cfg, out)

    click.echo("=== STEP 4: Repeat check (BLOSUM + is_repeat) ===")
    run_repeat_check(out)

    click.echo("=== STEP 5: Plotting report ===")
    plot_report(out)

    click.echo("=== STEP 6: Generating final report ===")
    generate_report(out, report_format)

    click.echo(f"\n🎉 FULL PIPELINE COMPLETE! 🎉")
    click.echo(f"Results in: {out}")
    click.echo("   • sample_qc_table.csv")
    click.echo("   • by_protein/*.csv.gz")
    click.echo("   • plots/")
    click.echo(f"   • NGS_Pipeline_Report.{report_format}")

if __name__ == "__main__":
    cli()
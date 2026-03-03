#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import jinja2
import numpy as np
from pathlib import Path
import sys
import argparse
import datetime
import json


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Interactive HTML Report Generator for OTU Pipeline."
    )
    parser.add_argument(
        "-f", "--filter_report",
        help="Path to Excel file 'report_filtering_...xlsx'",
        required=True, type=str
    )
    parser.add_argument(
        "-otu", "--otu_table",
        help="Path to Excel file 'Final_OTU_Table.xlsx'",
        required=True, type=str
    )
    parser.add_argument(
        "-o", "--output_html",
        help="Path for output HTML file (e.g., 'MyReport.html').",
        required=True, type=str
    )
    parser.add_argument(
        "-pn", "--project_name",
        help="[Optional] Project name to display on report.",
        type=str, default="Metagenomic Analysis Report"
    )
    return parser.parse_args()


def load_data(filter_path, otu_path):
    print(f"Loading filter data from: {filter_path}")
    try:
        df_filter = pd.read_excel(filter_path)
    except Exception as e:
        print(f"ERROR: Cannot read filter file: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"Loading OTU table from: {otu_path}")
    try:
        df_otu = pd.read_excel(otu_path, index_col=0)
    except Exception as e:
        print(f"ERROR: Cannot read OTU file: {e}", file=sys.stderr)
        sys.exit(1)
    return df_filter, df_otu


def get_kpis(df_filter, df_otu):
    print("Calculating KPIs...")
    total_seq = df_filter['total count'].sum()
    good_seq = df_filter['good count'].sum()
    if total_seq > 0:
        rejection_rate = f"{(total_seq - good_seq) / total_seq:.1%}"
    else:
        rejection_rate = "N/A"
    kpis = {
        'total_sequences': f"{total_seq:,}",
        'good_sequences': f"{good_seq:,}",
        'rejection_rate': rejection_rate,
        'otus_found': f"{len(df_otu):,}",
        'samples_processed': len(df_filter)
    }
    return kpis


def generate_filter_graphs(df_filter):
    print("Generating filter graphs...")
    df_filter_long = df_filter.melt(
        id_vars=['barcode'],
        value_vars=['good count', 'bad count'],
        var_name='Count Type',
        value_name='Sequences'
    )
    fig_counts = px.bar(
        df_filter_long, x='barcode', y='Sequences', color='Count Type',
        title="Filter Count Report (Good vs. Bad)",
        color_discrete_map={'good count': '#2ca02c', 'bad count': '#d62728'}
    )
    fig_counts_div = pio.to_html(fig_counts, full_html=False, include_plotlyjs='cdn')
    fig_lengths = px.bar(
        df_filter, x='barcode', y='mean length',
        title='Mean Sequence Length (post-filter)',
        labels={'mean length': 'Mean Length (bp)'}
    )
    fig_lengths.add_trace(go.Scatter(
        x=df_filter['barcode'], y=df_filter['median length'],
        mode='markers', name='Median Length',
        marker=dict(color='black', size=8)
    ))
    fig_lengths_div = pio.to_html(fig_lengths, full_html=False, include_plotlyjs=False)
    return fig_counts_div, fig_lengths_div


def get_positional_taxonomy(tax_str, position):
    try:
        if not isinstance(tax_str, str):
            return 'Unclassified'

        parts = tax_str.split(';')
        parts = [p for p in parts if p and not p.lower().startswith('unclassified')]

        if len(parts) >= abs(position):
            return parts[position]
        else:
            return 'Unclassified (Level not present)'

    except Exception:
        return 'Unclassified (Parsing Error)'


def aggregate_by_rank(df_otu, sample_cols, rank_name, position, top_n=15):
    print(f"Aggregating for {rank_name} (position {position})...")

    df_otu[rank_name] = df_otu['Taxonomy'].apply(
        lambda x: get_positional_taxonomy(x, position)
    )

    df_agg = df_otu.groupby(rank_name)[sample_cols].sum()

    df_agg['Total'] = df_agg.sum(axis=1)
    df_agg = df_agg.sort_values(by='Total', ascending=False)

    if len(df_agg) > top_n:
        df_top = df_agg.head(top_n - 1)
        df_others = df_agg.tail(len(df_agg) - (top_n - 1)).sum().to_frame().T
        df_others.index = [f'Others ({rank_name})']
        df_agg = pd.concat([df_top, df_others])

    df_agg = df_agg.drop(columns='Total')

    df_rel = df_agg.div(df_agg.sum(axis=0), axis=1) * 100

    df_rel.index.name = rank_name

    df_long = df_rel.reset_index().melt(
        id_vars=rank_name,
        var_name='Sample',
        value_name='Relative Abundance (%)'
    )
    df_long = df_long.rename(columns={rank_name: 'Taxonomy'})
    df_long['Rank'] = rank_name
    return df_long


def generate_otu_data_json(df_otu):
    print("Generating taxonomic data (Family, Genus, Species)...")

    if 'Taxonomy' not in df_otu.columns:
        print("WARNING: 'Taxonomy' column not found. Skipping OTU graphs.", file=sys.stderr)
        return "null", None

    sample_cols = [col for col in df_otu.columns if col != 'Taxonomy']

    df_family = aggregate_by_rank(df_otu, sample_cols, 'Family', -3, top_n=20)
    df_genus = aggregate_by_rank(df_otu, sample_cols, 'Genus', -2, top_n=15)
    df_species = aggregate_by_rank(df_otu, sample_cols, 'Species', -1, top_n=15)

    df_all_ranks = pd.concat([df_family, df_genus, df_species], ignore_index=True)

    composition_json_data = df_all_ranks.to_json(orient='records')

    print("Generating Heatmap graph...")
    df_counts_only = df_otu[sample_cols]
    top_25_otus = df_counts_only.sum(axis=1).nlargest(25).index
    df_top25 = df_counts_only.loc[top_25_otus]
    df_top25_log = np.log1p(df_top25)
    fig_heatmap = px.imshow(
        df_top25_log, title="Top 25 OTU Heatmap (Log-normalized Abundance)",
        labels=dict(x="Sample", y="OTU ID", color="Log(Count+1)"),
        x=df_top25_log.columns, y=df_top25_log.index,
        aspect="auto", height=700
    )
    fig_heatmap_div = pio.to_html(fig_heatmap, full_html=False, include_plotlyjs=False)

    return composition_json_data, fig_heatmap_div


def main():
    options = parse_arguments()

    filter_report_path = Path(options.filter_report)
    otu_table_path = Path(options.otu_table)
    output_html_path = Path(options.output_html)
    template_file = Path(__file__).parent / "template.html"

    if not template_file.exists():
        print("ERROR: 'template.html' file not found.", file=sys.stderr)
        sys.exit(1)

    df_filter, df_otu = load_data(filter_report_path, otu_table_path)
    kpis = get_kpis(df_filter, df_otu)
    graph_filter_counts_div, graph_filter_lengths_div = generate_filter_graphs(df_filter)
    composition_json_data, graph_otu_heatmap_div = generate_otu_data_json(df_otu)

    print("Converting tables to HTML...")
    table_filter_html = df_filter.to_html(
        classes="table table-striped table-hover table-sm",
        index=False, border=0
    )
    table_otu_html = df_otu.to_html(
        classes="table table-striped table-hover table-sm",
        index=True, border=0, table_id="otu-table"
    )

    print("Rendering final HTML report...")
    template_loader = jinja2.FileSystemLoader(searchpath=Path(__file__).parent)
    template_env = jinja2.Environment(loader=template_loader)
    template = template_env.get_template(template_file.name)

    context = {
        "project_name": options.project_name,
        "report_date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "kpis": kpis,
        "graph_filter_counts_div": graph_filter_counts_div,
        "graph_filter_lengths_div": graph_filter_lengths_div,
        "table_filter_html": table_filter_html,
        "composition_json_data": composition_json_data,
        "graph_otu_heatmap_div": graph_otu_heatmap_div,
        "table_otu_html": table_otu_html
    }

    try:
        html_content = template.render(context)
        with open(output_html_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        print(f"\n--- HTML Report Generated Successfully! ---")
        print(f"Open this file in your browser: {output_html_path.resolve()}")

    except Exception as e:
        print(f"ERROR during template rendering: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

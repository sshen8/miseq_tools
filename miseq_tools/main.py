import argparse
from .samplesheet import format_samplesheet
from .quant_kapa import kapaquant
from .quant_qubit import qubitquant
from .pooling import pooling
from .quant_combine import quant_combine
from .demux_stats import demux
import logging
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Tools for running Miseq", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # subparser for sample sheet
    subparsers = parser.add_subparsers(required=True)
    parser.add_argument("--mpl-style", help="Matplotlib style to use", default="simon")
    parser.add_argument("--log-level", help="Log level", default="INFO", choices=logging._nameToLevel.keys())
    parser_samplesheet = subparsers.add_parser("sheet", help="Format sample sheet for Miseq")
    parser_samplesheet.set_defaults(func=format_samplesheet)
    parser_samplesheet.add_argument("fname_in", help="Input file")
    parser_samplesheet.add_argument("-o", help="Output file", dest="fname_out", default="samplesheet.csv")

    parser_kapa = subparsers.add_parser("kapa", help="Analyze qPCR library quantification data")
    parser_kapa.set_defaults(func=kapaquant)
    parser_kapa.add_argument("kapafolder", help="Folder containing qPCR data")
    parser_kapa.add_argument("samplesheet", help="Sample sheet to use")
    parser_kapa.add_argument("--dilution", help="Dilution factor of samples", type=float, default=1e4)
    parser_kapa.add_argument("--standard-bp", help="Amplicon size (bp) of standards. 452 for KAPA, 399 for NEB.", type=int, default=399)

    parser_qubit = subparsers.add_parser("qubit", help="Analyze Qubit quantification data")
    parser_qubit.add_argument("samplesheet", help="Sample sheet to use")
    parser_qubit.set_defaults(func=qubitquant)

    parser_combine = subparsers.add_parser("combine", help="Combine KAPA and Qubit quantification data")
    parser_combine.add_argument("--kapa", help="KAPA quantification data", dest="kapa_fname", default="quant_kapa.csv")
    parser_combine.add_argument("--qubit", help="Qubit quantification data", dest="qubit_fname", default="quant_qubit.csv")
    parser_combine.set_defaults(func=quant_combine)

    parser_pooling = subparsers.add_parser("pool", help="Figure out pooling")
    parser_pooling.add_argument("samplesheet", help="Sample sheet to use")
    parser_pooling.add_argument("quant_csv", help="Quantification data")
    parser_pooling.add_argument("--pools", help="Number of pools", type=int, default=3)
    parser_pooling.set_defaults(func=pooling)

    parser_pre = subparsers.add_parser("pre", help="Full pre-Miseq pipeline: includes sheet, kapa, qubit, combine, and pool")
    parser_pre.add_argument("samplesheet", help="Sample sheet to use")
    parser_pre.add_argument("kapafolder", help="Folder containing KAPA data")
    def pipeline_pre(**kwargs):
        format_samplesheet(fname_in=kwargs['samplesheet'], fname_out='samplesheet.csv')
        kapaquant(kapafolder=kwargs['kapafolder'], samplesheet=kwargs['samplesheet'], dilution=1e4, standard_bp=399)
        qubitquant(samplesheet=kwargs['samplesheet'])
        quant_combine(kapa_fname="quant_kapa.csv", qubit_fname="quant_qubit.csv")
        pooling(samplesheet=kwargs['samplesheet'], quant_csv="quant_combined.csv", pools=3)
    parser_pre.set_defaults(func=pipeline_pre)

    parser_demux = subparsers.add_parser("demux", help="Demuxing stats")
    parser_demux.set_defaults(func=demux)
    parser_demux.add_argument("samplesheet", help="Sample sheet to use")
    parser_demux.add_argument("stats", help="Stats.json file from Miseq")

    args = parser.parse_args()

    plt.style.use(args.mpl_style)
    del args.mpl_style
    logging.basicConfig(format='%(levelname)-10s%(message)s', level=args.log_level)
    del args.log_level
    kwargs = vars(args)

    args.func(**{k: v for k, v in kwargs.items() if k != "func"})

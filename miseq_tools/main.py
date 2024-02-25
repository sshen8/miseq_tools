import argparse
from .samplesheet import format_samplesheet
from .quant_kapa import kapaquant
from .quant_qubit import qubitquant
from .pooling import pooling
from .quant_combine import quant_combine
import logging

def main():
    logging.basicConfig(format='%(levelname)-10s%(message)s', level=logging.INFO)

    parser = argparse.ArgumentParser(description="Tools for running Miseq")
    # subparser for sample sheet
    subparsers = parser.add_subparsers()
    parser_samplesheet = subparsers.add_parser("sheet", help="Format sample sheet for Miseq")
    parser_samplesheet.set_defaults(func=format_samplesheet)
    parser_samplesheet.add_argument("fname_in", help="Input file")
    parser_samplesheet.add_argument("-o", help="Output file", dest="fname_out", default="samplesheet.csv")

    parser_kapa = subparsers.add_parser("kapa", help="Analyze KAPA quantification data")
    parser_kapa.set_defaults(func=kapaquant)
    parser_kapa.add_argument("kapafolder", help="Folder containing KAPA data")
    parser_kapa.add_argument("samplesheet", help="Sample sheet to use")
    parser_kapa.add_argument("--dilution", help="Dilution factor for KAPA quantification", type=float, default=1e4)

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
    parser_pooling.set_defaults(func=pooling)

    parser_pre = subparsers.add_parser("pre", help="Full pre-Miseq pipeline: includes sheet, kapa, qubit, combine, and pool")
    parser_pre.add_argument("samplesheet", help="Sample sheet to use")
    parser_pre.add_argument("kapafolder", help="Folder containing KAPA data")
    def pipeline_pre(**kwargs):
        format_samplesheet(fname_in=kwargs['samplesheet'], fname_out='samplesheet.csv')
        kapaquant(kapafolder=kwargs['kapafolder'], samplesheet=kwargs['samplesheet'], dilution=1e4)
        qubitquant(samplesheet=kwargs['samplesheet'])
        quant_combine(kapa_fname="quant_kapa.csv", qubit_fname="quant_qubit.csv")
        pooling(samplesheet=kwargs['samplesheet'], quant_csv="quant_combined.csv")
    parser_pre.set_defaults(func=pipeline_pre)

    args = parser.parse_args()
    kwargs = vars(args)
    args.func(**{k: v for k, v in kwargs.items() if k != "func"})

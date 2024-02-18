import sys
import os
from pathlib import Path
from dp.launching.typing import BaseModel, Field, OutputDirectory,InputFilePath, Int,  Union, String, Literal, Field, Enum
from dp.launching.cli import SubParser,default_minimal_exception_handler,run_sp_and_exit,to_runner
import pandas as pd
def retrosynthesis_func(mol_file, id_tag, output_dir):
    os.system("mkdir -p %s" % (output_dir))

    os.system( f" export PYTHONPATH=\"$PYTHONPATH:/root/ASKCOS\" && /root/miniconda3/envs/ASKCOS_0.3.1/bin/python \
    /root/ASKCOS/makeit/retrosynthetic/tree_builder_bohrium_task.py  --input_csv_file {mol_file}  \
    --output_dir {output_dir}  --ID_tag {id_tag} ")




class Options(BaseModel):
    mol_file: InputFilePath = Field(description="""CSV ,the CSV file has 2 columns, "smiles" or "SMILES" and id tag which will be indicated on next step .""")
    id_tag:String=Field(description="The name of the ID tag in input CSV file, which could not have space character")
    output_dir: OutputDirectory = Field(
        default="./outputs"
    )  


class global_opt(Options, BaseModel):
    ...


def runner(opts:  global_opt) -> int:
    retrosynthesis_func(mol_file=opts.mol_file, id_tag=opts.id_tag, output_dir=opts.output_dir)


def to_parser():
    return to_runner(
        global_opt,
        runner,
        version="0.1.0",
        exception_handler=default_minimal_exception_handler,
    )


if __name__ == '__main__':
    import sys
    to_parser()(sys.argv[1:])

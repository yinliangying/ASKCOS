import sys
import os
from pathlib import Path
from dp.launching.typing import BaseModel, Field,Boolean, OutputDirectory,InputFilePath, Int,  Union, String, Literal, Field, Enum
from dp.launching.cli import SubParser,default_minimal_exception_handler,run_sp_and_exit,to_runner

import pandas as pd
def retrosynthesis_func(mol_file, predicting_reaction_condition, output_dir):
    os.system("mkdir -p %s" % (output_dir))

    os.system( f" export PYTHONPATH=$PYTHONPATH:/root/Uni-Electrolyte/retrosynthesis/ASKCOS2:/root/Uni-Electrolyte/retrosynthesis/ASKCOS2/askcos "
               f"&& cd /root/Uni-Electrolyte/retrosynthesis/g2gretro/ && nohup sh -x start_server.sh& && sleep 10"
               f"&& /root/miniconda3/envs/ASKCOS_0.3.1/bin/python /root/Uni-Electrolyte/retrosynthesis/ASKCOS2/tree_builder_bohrium_task.py  "
               f"--input_csv_file {mol_file}  --output_dir {output_dir}  "
               f"--predicting_reaction_condition {predicting_reaction_condition} "
               )

    if not os.path.exists(f"{output_dir}/finish"):
        raise Exception("task failed")


class Options(BaseModel):
    mol_file: InputFilePath = Field(description="""CSV ,the CSV file has 2 columns, "smiles" or "SMILES" .""")
    predicting_reaction_condition:Boolean
    output_dir: OutputDirectory = Field(
        default="./outputs"
    )



class global_opt(Options, BaseModel):
    ...


def runner(opts:  global_opt) -> int:
    retrosynthesis_func(mol_file=opts.mol_file, predicting_reaction_condition=opts.predicting_reaction_condition, output_dir=opts.output_dir)


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

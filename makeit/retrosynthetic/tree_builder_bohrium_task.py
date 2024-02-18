import makeit.global_config as gc
from makeit.retrosynthetic.tree_builder import TreeBuilder
from argparse import ArgumentParser
import os
import pandas as pd
import json
import time
from collections import deque
import rdkit
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from PIL import Image

parser = ArgumentParser()
parser.add_argument('--input_csv_file', type=str)
parser.add_argument('--output_dir', type=str)
parser.add_argument("--ID_tag", type=str)
args = parser.parse_args()

input_df = pd.read_csv(args.input_csv_file)
ID_tag=args.ID_tag
output_dir = args.output_dir
# if os.path.exists(output_dir):
#     os.system("rm -rf %s" % (output_dir))
os.system("mkdir %s" % (output_dir))


celery = False
st = time.time()
treeBuilder = TreeBuilder(celery=celery, mincount=25, mincount_chiral=10)
print("treeBuilder spend", time.time() - st)

result_dict = {}
os.system("mkdir %s/tmp"%(output_dir))
for k, v in input_df.iterrows():
    try:
        smiles = v["smiles"]
    except:
        smiles=v["SMILES"]
    ID = v[ID_tag]
    result_dict[ID] = {"smiles": smiles}
    status, paths = treeBuilder.get_buyable_paths(smiles, max_depth=4, template_prioritization=gc.relevance,
                                                  precursor_prioritization=gc.relevanceheuristic, nproc=2,
                                                  expansion_time=60, max_trees=5, max_ppg=10,
                                                  max_branching=25, apply_fast_filter=True, filter_threshold=0.75,
                                                  min_chemical_history_dict={'as_reactant': 5, 'as_product': 1,
                                                                             'logic': 'none'})
    result_dict[ID]["paths"] = paths
    result_dict[ID]["status"] = status
    # if len(paths) < 3:
    #     result_dict[ID]["paths"] = paths
    # else:
    #     result_dict[ID]["paths"] = paths[:3]

output_json_fp=open("%s/output.json"%(output_dir),"w")
json.dump(result_dict, output_json_fp)

# import os
# import pandas as pd
# import json
# import time
# from collections import deque
# import rdkit
# from rdkit.Chem import AllChem
# from rdkit.Chem import Draw
# from PIL import Image

# tmp_dir="./tmp_dir"
# #tmp_dir="./data/rxn_draw"
# if os.path.exists(tmp_dir):
#     os.system("rm -rf %s"%(tmp_dir))
# os.system("mkdir %s"%(tmp_dir))
# result_dict=json.load(open("./output_json_file.json"))


q = deque()
for molID in result_dict:
    os.system("mkdir %s/molecule_%s"%(output_dir,molID))
    for path_id, path_dict in enumerate(result_dict[molID]["paths"]): #广度遍历tree

        head = path_dict
        q.append(head)
        rxn_smiles_list = []
        while (len(q) != 0):
            point = q.popleft()
            if ">" in point["smiles"]:
                rxn_smiles_list.append(point["smiles"])
            for cp in point["children"]:
                q.append(cp)
        # print()

        if len(rxn_smiles_list) == 0:
            continue
        elif len(rxn_smiles_list) == 1:
            # continue
            rxn = AllChem.ReactionFromSmarts(rxn_smiles_list[0], useSmiles=True)
            img = Draw.ReactionToImage(rxn, subImgSize=(800, 300))
            img.save('%s/molecule_%s/pathway_%s_%s.png' % (output_dir,molID, molID, path_id))
            # img.show()
            # exit()
            # d2d = Draw.MolDraw2DCairo(800, 300)
            # d2d.DrawReaction(rxn)
            # png = d2d.GetDrawingText()
            # open('%s/result_%s_%s.jpg'%(tmp_dir,molID,path_id), 'wb+').write(png)
        else:
            pil_img_list = []
            width = 800
            height = 300
            width_mol = 200
            for rxn_idx, rxn_smiles in enumerate(reversed(rxn_smiles_list)):
                rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
                img = Draw.ReactionToImage(rxn, subImgSize=(width_mol, height))  # 每个分子的尺寸 反应箭头也按一个分子算
                img = img.resize((width, height))
                img.save("%s/tmp/tmp_%s_%s_%s.png" % (output_dir, molID, path_id, rxn_idx))
                # d2d = Draw.MolDraw2DCairo(800, 300)
                # d2d.DrawReaction(rxn)
                # png = d2d.GetDrawingText()
                # open("%s/tmp_%s_%s_%s.png"%(tmp_dir,molID,path_id,rxn_idx), 'wb+').write(png)
                # img=Draw.ReactionToImage(rxn,subImgSize=(200, 200))
                # img.save("./data/rxn_draw/%s_%s_%s.png"%(molID,path_id,rxn_idx))

                pil_img_list.append(img)

            # 创建一个空白画布，用于拼接图片
            result_width = width  # 图片拼接在一起
            result_height = height * len(pil_img_list)
            result = Image.new(pil_img_list[0].mode, (result_width, result_height))  #

            # 在画布上拼接图片
            for img_idx, img in enumerate(pil_img_list):
                result.paste(img, (0, height * img_idx))

            # 保存拼接后的图片
            result.save('%s/molecule_%s/pathway_%s_%s.png' % (output_dir,molID, molID, path_id))

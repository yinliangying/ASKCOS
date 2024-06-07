import sys

import makeit.global_config as gc
from makeit.retrosynthetic.tree_builder import TreeBuilder
from makeit.synthetic.context.neuralnetwork import  NeuralNetContextRecommender
from argparse import ArgumentParser
import os
import pandas as pd
import json
import traceback
import time
from collections import deque
import rdkit
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont

cont = NeuralNetContextRecommender()
cont.load_nn_model(model_path=gc.NEURALNET_CONTEXT_REC['model_path'],
                   info_path=gc.NEURALNET_CONTEXT_REC['info_path'],
                   weights_path=gc.NEURALNET_CONTEXT_REC['weights_path'])
width = 1200
height_mol = 300
width_mol = 300

font_file="DejaVuSans.ttf" # "arial.ttf" for win
def output_condiction_picture(rxn_smiles,condition_result):

    pil_img_list = []


    # 反应的picture
    rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    img = Draw.ReactionToImage(rxn, subImgSize=(width_mol, height_mol))
    img = img.resize((width, height_mol))
    pil_img_list.append(img)

    #拼接条件title
    condition_part_img_list = []
    for sth in ("Temperature","Solvent", "Reagents", "Catalyst"):
        sth_img = Image.new(pil_img_list[0].mode, (width_mol, height_mol), color=(0, 0, 0))  #
        draw = ImageDraw.Draw(sth_img)
        draw.rectangle([(0, 0), sth_img.size], fill=(255, 255, 255))
        text = sth
        position = (50,250)  # 文字的起始位置 (x, y)
        font = ImageFont.truetype(font_file, 24)  # 使用指定字体和大小
        color = (0, 0, 0)  # 文字颜色，RGB 格式
        draw.text(position, text, font=font, fill=color)
        condition_part_img_list.append(sth_img)
    # 创建一个空白画布，用于拼接condition图片
    condition_img = Image.new(pil_img_list[0].mode, (width, height_mol), color=(0, 0, 0))  #
    # 在画布上拼接图片
    for img_idx, img in enumerate(condition_part_img_list):
        condition_img.paste(img, (width_mol * img_idx, 0))
    pil_img_list.append(condition_img)


    # condition_result=json.loads("""
    # [[102.30387878417969, "C1COCCO1", "CCN(CC)CC", "Reaxys Name (1,1'-bis(diphenylphosphino)ferrocene)palladium(II) dichloride", NaN, NaN], [104.92787170410156, "C1COCCO1", "CCN(CC)CC", "Cl[Pd](Cl)([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1", NaN, NaN], [99.1409912109375, "Cc1ccccc1", "CCN(CC)CC", "Cl[Pd](Cl)([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1", NaN, NaN], [76.38555908203125, "C1CCOC1", "CCN(CC)CC", "Cl[Pd](Cl)([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1", NaN, NaN], [95.92562103271484, "Cc1ccccc1", "CCN(CC)CC", "Reaxys Name (1,1'-bis(diphenylphosphino)ferrocene)palladium(II) dichloride", NaN, NaN], [75.68882751464844, "C1CCOC1", "CCN(CC)CC", "Reaxys Name (1,1'-bis(diphenylphosphino)ferrocene)palladium(II) dichloride", NaN, NaN], [93.39191436767578, "C1COCCO1", "", "Reaxys Name (1,1'-bis(diphenylphosphino)ferrocene)palladium(II) dichloride", NaN, NaN], [97.8741226196289, "C1COCCO1", "CC(=O)[O-].[K+]", "Reaxys Name (1,1'-bis(diphenylphosphino)ferrocene)palladium(II) dichloride", NaN, NaN], [95.84452819824219, "C1COCCO1", "[MgH2]", "Cl[Pd](Cl)([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1", NaN, NaN], [67.86063385009766, "C1CCOC1", "[MgH2]", "Cl[Pd](Cl)([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1", NaN, NaN]]""")
    for condition_list in condition_result:
        temperature, solvent, reagents, catalyst, _, _ = condition_list
        condition_part_img_list = []

        temperature_img = Image.new(pil_img_list[0].mode, (width_mol, height_mol), color=(255, 255, 255))  #
        draw = ImageDraw.Draw(temperature_img)
        draw.rectangle([(0, 0), temperature_img.size], fill=(255, 255, 255))
        text = "%s°C" % (int(temperature))
        position = (0,100)  # 文字的起始位置 (x, y)
        font = ImageFont.truetype(font_file, 24)  # 使用指定字体和大小
        color = (0, 0, 0)  # 文字颜色，RGB 格式
        draw.text(position, text, font=font, fill=color)
        condition_part_img_list.append(temperature_img)

        for sth in (solvent, reagents, catalyst):
            mol = AllChem.MolFromSmiles(sth)
            if not mol:
                sth_img = Image.new(pil_img_list[0].mode, (width_mol, height_mol), color=(255, 255, 255))  #
                draw = ImageDraw.Draw(sth_img)
                draw.rectangle([(0, 0), sth_img.size], fill=(255,255,255))
                text=""
                for char_idx ,char in enumerate(sth):
                    if char_idx%20==0:
                        text+="\n"
                    text+=char
                position = (0,0)  # 文字的起始位置 (x, y)
                font = ImageFont.truetype(font_file, 24)  # 使用指定字体和大小
                color = (0,0,0)  # 文字颜色，RGB 格式
                draw.text(position, text, font=font, fill=color)
            else:
                sth_img = Draw.MolToImage(mol)
            condition_part_img_list.append(sth_img)

        # 创建一个空白画布，用于拼接condition图片
        condition_img = Image.new(pil_img_list[0].mode, (width, height_mol), color=(255, 255, 255))  #
        # 在画布上拼接图片
        for img_idx, img in enumerate(condition_part_img_list):
            condition_img.paste(img, (width_mol * img_idx, 0))

        pil_img_list.append(condition_img)

    # 创建一个空白画布，用于拼接图片
    result_width = width  # 图片拼接在一起
    result_height = height_mol * len(pil_img_list)
    condition_img_result = Image.new(pil_img_list[0].mode, (result_width, result_height),color=(255, 255, 255))  #
    # 在画布上拼接图片
    for img_idx, img in enumerate(pil_img_list):
        condition_img_result.paste(img, (0, height_mol * img_idx))
    #result.show()
    return condition_img_result



def draw_pathway(result_dict,args,output_dir):


    for mol_idx in result_dict:
        # if len(result_dict[mol_idx]["paths"])!=0:
        #     os.system("mkdir %s/molecule_%s" % (output_dir, mol_idx))

        for path_id, path_dict in enumerate(result_dict[mol_idx]["paths"]):  # 广度遍历tree
            q = deque()
            head = path_dict
            q.append(head)
            rxn_info_list = []
            price_dict = {}
            while (len(q) != 0):
                point = q.popleft()
                if ">" in point["smiles"]:

                    #rxn_condition_dir_name="molecule_%s_pathway_%s_rxn_id_%s_condition" % (mol_idx, path_id, len(rxn_info_list))
                    rxn_info_list.append((point["smiles"],"",point["condition_result"]))
                if "ppg" in point:
                    price_dict[point["smiles"]] = point["ppg"]
                for cp in point["children"]:
                    q.append(cp)
            # print()

            if len(rxn_info_list) == 0:
                continue

            else:
                pil_img_list = []

                if True:#args.predicting_reaction_condition:
                    os.system(
                        "mkdir -p %s/molecule_%s/pathway_%s_%s_condition" % (output_dir, mol_idx, mol_idx, path_id))
                for rxn_idx, rxn_info in enumerate(reversed(rxn_info_list)):
                    rxn_smiles=rxn_info[0]
                    condition_result=rxn_info[2]
                    rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
                    img = Draw.ReactionToImage(rxn, subImgSize=(width_mol, height_mol))  # 每个分子的尺寸 反应箭头也按一个分子算
                    img = img.resize((width, height_mol))

                    reactant_smiles_str=rxn_smiles.split(">")[0]
                    reactant_smiles_list=reactant_smiles_str.split(".")
                    tmp_width_mol=width/(2+len(reactant_smiles_list))
                    for reactant_idx,reactant_smiles in enumerate(reactant_smiles_list):
                        if reactant_smiles in price_dict:
                            price = float(price_dict[reactant_smiles])
                            if price==0:
                                price+=0.1
                            text="%.1f/g"%(price)
                            draw = ImageDraw.Draw(img)

                            position = (100+reactant_idx*tmp_width_mol,height_mol-50 )  # 文字的起始位置 (x, y)
                            font = ImageFont.truetype(font_file, 24)  # 使用指定字体和大小
                            color = (0, 0, 0)  # 文字颜色，RGB 格式
                            draw.text(position, text, font=font, fill=color)

                    img.save("%s/tmp/tmp_%s_%s_%s.png" % (output_dir, mol_idx, path_id, rxn_idx))
                    # d2d = Draw.MolDraw2DCairo(800, 300)
                    # d2d.DrawReaction(rxn)
                    # png = d2d.GetDrawingText()
                    # open("%s/tmp_%s_%s_%s.png"%(tmp_dir,molID,path_id,rxn_idx), 'wb+').write(png)
                    # img=Draw.ReactionToImage(rxn,subImgSize=(200, 200))
                    # img.save("./data/rxn_draw/%s_%s_%s.png"%(molID,path_id,rxn_idx))

                    pil_img_list.append(img)

                    # condition prediction
                    if True:#args.predicting_reaction_condition:
                        condition_img  = output_condiction_picture(rxn_smiles,condition_result)
                        condition_img.save(
                            "%s/molecule_%s/pathway_%s_%s_condition/rxn_%s_condition.png" % (
                            output_dir, mol_idx, mol_idx, path_id, rxn_idx))

                # 创建一个空白画布，用于拼接图片
                result_width = width  # 图片拼接在一起
                result_height = height_mol * len(pil_img_list)
                result = Image.new(pil_img_list[0].mode, (result_width, result_height), color=(255, 255, 255))  #

                # 在画布上拼接图片
                for img_idx, img in enumerate(pil_img_list):
                    result.paste(img, (0, height_mol * img_idx))

                # 保存拼接后的图片
                result.save('%s/molecule_%s/pathway_%s_%s.png' % (output_dir, mol_idx, mol_idx, path_id))




def main():
    parser = ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--output_dir', type=str)
    parser.add_argument("--predicting_reaction_condition", type=str)
    args = parser.parse_args()
    if args.predicting_reaction_condition == "True":
        args.predicting_reaction_condition = True
    else:
        args.predicting_reaction_condition = False

    # ID_tag=args.ID_tag
    output_dir = args.output_dir
    if os.path.exists(output_dir):
        os.system("rm -rf %s" % (output_dir))
    os.system("mkdir %s" % (output_dir))

    celery = False
    st = time.time()
    treeBuilder = TreeBuilder(celery=celery, mincount=25, mincount_chiral=10)
    print("treeBuilder spend", time.time() - st)



    result_dict = {}
    os.system("mkdir %s/tmp" % (output_dir))
    for idx, smiles in enumerate(open(args.input_file)):

        #
        # import pdb
        # pdb.set_trace()
        # smiles="CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2"
        smiles = smiles.strip()
        result_dict[idx] = {"smiles": smiles}
        try:
            status, paths = treeBuilder.get_buyable_paths(smiles, max_depth=4, template_prioritization=gc.relevance,
                                                          precursor_prioritization=gc.relevanceheuristic, nproc=2,
                                                          expansion_time=60, max_trees=15, max_ppg=10,
                                                          max_branching=25, apply_fast_filter=True,
                                                          filter_threshold=0.75,
                                                          min_chemical_history_dict={'as_reactant': 5, 'as_product': 1,
                                                                                     'logic': 'none'})
        except:
            print(smiles, file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            continue
        print("%s %s complete" % (idx, smiles), file=sys.stderr)

        result_dict[idx]["paths"] = paths
        result_dict[idx]["status"] = status
        # if len(paths) < 3:
        #     result_dict[ID]["paths"] = paths
        # else:
        #     result_dict[ID]["paths"] = paths[:3]

        for path_id, path_dict in enumerate(result_dict[idx]["paths"]):
            # 广度遍历tree
            q = deque()
            q.append(path_dict)
            while (len(q) != 0):
                point = q.popleft()
                if ">" in point["smiles"]:
                    rxn_smiles=point["smiles"]
                    condition_result = cont.get_n_conditions(rxn_smiles, 10, with_smiles=False)
                    if condition_result:
                        point["condition_result"] = condition_result
                for cp in point["children"]:
                    q.append(cp)


    output_json_fp = open("%s/output.json" % (output_dir), "w")
    json.dump(result_dict, output_json_fp)

    draw_pathway(result_dict,args,output_dir)


    with open("%s/finish" % (output_dir), "w") as fp:
        print("OK", file=fp)


if __name__=="__main__":
    main()

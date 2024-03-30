import http.client
from urllib.parse import quote
import json
# 发送 GET 请求到服务器
def send_get_request(path, query_params=None):
    conn = http.client.HTTPConnection("localhost", 8000)
    url = path
    if query_params:
        url += "?" + "&".join(["%s=%s"%(key,value) for key, value in query_params.items()])
    conn.request("GET", url)
    response = conn.getresponse()
    print("Response Status: %s"%(response.status))
    return response.read().decode('utf-8')

def get_ASKCOS_one_step_retro_topN(smiles,topN):
    json_str = json.dumps(
        {"smiles_list": [smiles]})
    encoded_value = quote(json_str)
    json_str=send_get_request("/", {"encoded_value": encoded_value})
    result_list=json.loads(json_str)
    condicate_smiles_list=result_list[0]

    return_list=[]
    for idx,r_smiles in enumerate(condicate_smiles_list):
        r_smiles_list=r_smiles.split(".")
        data={'necessary_reagent':"",
              'tforms':[],
              'plausibility':0.9,
              'template_score':0.01,
              'rank':idx+1,
              'num_examples':100,
              'smiles_split':r_smiles_list,
         'smiles':r_smiles, 'score':-1000*idx}
        return_list.append(data)
        if len(return_list)==topN:
            break
    return return_list


if __name__=="__main__":
    aa=get_ASKCOS_one_step_retro_topN("CC1(C)CC2COC1C2")
    print(aa)
    exit()

    # 测试 GET 请求


    print("Sending GET request...")
    json_str=json.dumps({"smiles_list": ["C=C1COC1CC=CC","Cc1ccc2oc-2co1","C=C1C#CC(CC)=CO1","CCC1COCOC1C","C=C1CCOC(C)O1",
                                           "CCOC=COCOCC","CC1(C)CC2COC1C2","C=CCC(=C)OC(C)OC","CCOC(COC)OCC","C=CC1=CCC#COC1",
                                           "C=COC=CC(=C)COC","CC12CC13COCC3O2","CCC(C)=CCC(C)=O","COC12C#CCCC1C2","C#CC1(C)CC#CC1=O"]})
    encoded_value = quote(json_str)
    print(send_get_request("/", {"encoded_value":encoded_value}))
    print()
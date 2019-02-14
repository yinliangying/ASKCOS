from django.http import JsonResponse
from makeit import global_config as gc
from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors as get_top_precursors_c

def singlestep(request):
    resp = {}
    target = request.GET.get('target')
    template_prioritization = request.GET.get('template_prioritization', gc.relevance)
    precursor_prioritization = request.GET.get('precursor_prioritization', gc.relevanceheuristic)
    mincount = request.GET.get('mincount', 0)
    num_templates = int(request.GET.get('num_templates', 100))
    apply_fast_filter = request.GET.get('apply_fast_filter', True)
    filter_threshold = float(request.GET.get('filter_threshold', 0.75))
    res = get_top_precursors_c.delay(target, template_prioritization, precursor_prioritization, mincount=mincount, apply_fast_filter=apply_fast_filter, filter_threshold=filter_threshold)
    (smiles, precursors) = res.get(100)
    resp['precursors'] = precursors
    return JsonResponse(resp)
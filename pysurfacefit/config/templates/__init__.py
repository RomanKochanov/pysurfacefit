from collections import OrderedDict

from ..base import Config
from . import general, data, model, fit, stat, codegen, plotting, multicut, calc

template_modules_list = general, data, model, fit, stat, codegen, plotting, multicut, calc

template_modules_dict = OrderedDict()
for mod in template_modules_list:
    template_modules_dict[mod.__name__.split('.')[-1]] = mod

config = Config(*[getattr(module,'Default') for module in template_modules_list])

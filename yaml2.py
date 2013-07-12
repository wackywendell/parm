"""Make PyYAML output an OrderedDict.

It will do so fine if you use yaml.dump(), but that generates ugly, 
non-standard YAML code.

To use yaml.safe_dump(), you need the following.
"""

from collections import OrderedDict
import yaml
from yaml import *

def represent_odict(dump, tag, mapping, flow_style=None):
    """Like BaseRepresenter.represent_mapping, but does not issue the sort().
    """
    value = []
    node = yaml.MappingNode(tag, value, flow_style=flow_style)
    if dump.alias_key is not None:
        dump.represented_objects[dump.alias_key] = node
    best_style = True
    if hasattr(mapping, 'items'):
        mapping = mapping.items()
    for item_key, item_value in mapping:
        node_key = dump.represent_data(item_key)
        node_value = dump.represent_data(item_value)
        if not (isinstance(node_key, yaml.ScalarNode) and not node_key.style):
            best_style = False
        if not (isinstance(node_value, yaml.ScalarNode) and not node_value.style):
            best_style = False
        value.append((node_key, node_value))
    if flow_style is None:
        if dump.default_flow_style is not None:
            node.flow_style = dump.default_flow_style
        else:
            node.flow_style = best_style
    return node

yaml.SafeDumper.add_representer(OrderedDict,
    lambda dumper, value: represent_odict(dumper, u'tag:yaml.org,2002:map', value))


# see http://hooke.googlecode.com/svn-history/r291/branches/wtk/hooke/util/yaml.py
import numpy

def int_representer(dumper, data):
    return dumper.represent_int(data)
yaml.SafeDumper.add_representer(numpy.int32, int_representer)
yaml.SafeDumper.add_representer(numpy.dtype(numpy.int32), int_representer)

def float_representer(dumper, data):
    return dumper.represent_float(float(data))
yaml.SafeDumper.add_representer(numpy.float32, float_representer)
yaml.SafeDumper.add_representer(numpy.float64, float_representer)


def array_representer(dumper, data, flow_style=None):
    return dumper.represent_sequence(u'tag:yaml.org,2002:seq', list(data), flow_style=flow_style)
yaml.SafeDumper.add_representer(numpy.ndarray, array_representer)

#~ def represent_np_float(dump, tag, item, flow_style=None):
    #~ return dump.represent_data(float(item))

"""Make PyYAML load a file that uses lists as mapping keys.
 
You would otherwise get an error like:
 
yaml.constructor.ConstructorError: while constructing a mapping
in "eiscp-commands.yaml", line 95, column 7
found unacceptable key (unhashable type: 'list')
"""
 
def construct_tuple(loader, node):
    return tuple(yaml.SafeLoader.construct_sequence(loader, node))

yaml.SafeLoader.add_constructor(u'tag:yaml.org,2002:seq', construct_tuple)

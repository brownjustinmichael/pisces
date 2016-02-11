import os
import yaml
import collections

def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

def process(value, dictionary=None, current_key=""):
    """
    Process a YAML nested dictionary into a flat dictionary with different levels becoming '__'-separated keys
    """
    if dictionary is None:
        dictionary = {}
    if isinstance(value, dict):
        for key in value:
            if current_key:
                composite_key="__".join([current_key, key])
            else:
                composite_key=key
            process(value[key], dictionary, composite_key)
    else:
        dictionary[current_key]=value
    return dictionary

def unprocess(value, dictionary=None):
    if dictionary is None:
        dictionary = {}
    for key in value:
        keys = key.split("__")
        dict = dictionary
        while len(keys) > 1:
            if keys[0] not in dict:
                dict[keys[0]] = {}
            dict = dict[keys[0]]
            keys.pop(0)
        dict[keys[0]] = value[key]
    return dictionary

def Configuration(file=None, default=None, **kwargs):
    """
    Returns a dictionary setup object for use with code objects

    :param file: The file from which the configuration should be loaded
    :param kwargs: Any additional parameters that should be added can be added as kwargs
    :return: The configuration dictionary
    :rtype: `dict`
    """
    # Load the default configuration
    if default is not None:
        tmp = yaml.load(open(default))
    else:
        tmp = yaml.load(open(os.path.join(os.path.dirname(__file__), "../src/defaults.yaml")))

    # Load any additional keys from the given file
    try:
        file = open(file)
    except Exception as e:
        pass
    if file is not None:
        tmp = update(tmp, yaml.load(file))

    # Update any additional arguments provided to the function
    tmp.update (kwargs)

    # Add some defaults if they haven't been specified
    if "np" not in tmp:
        tmp ["np"] = 1
    if "wd" not in tmp:
        tmp ["wd"] = os.getcwd()

    return tmp
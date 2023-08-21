import yaml
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '')] + sys.path

with open('version.h', 'r') as f:
    for line in f:
        if line.startswith('#define'):
            version = line.split()[2].strip('"')


data = {
    "package": {
        "name": "kma",
        "version": version
    },
    "source": {
        "url": "https://bitbucket.org/genomicepidemiology/kma/get/{}.tar.gz".format(version),
    },
    "build": {
        "number": 0
    },
    "requirements": {
        "build": [
            "make",
            "{{ compiler('c') }}"
        ],
        "host": [
            "zlib"
        ],
        "run": [
            "zlib"
        ]
    },
    "about": {
        "home": "https://bitbucket.org/genomicepidemiology/kma",
        "summary": "KMA is mapping a method designed to map raw reads directly against redundant databases, in an ultra-fast manner using seed and extend.",
        "license": "Apache-2.0"
    },
    "extra": {
        "identifiers": {
            "doi": "10.1186/s12859-018-2336-6",
        }
    }
}

# Convert the data to YAML and print it
yaml_str = yaml.dump(data, sort_keys=False)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
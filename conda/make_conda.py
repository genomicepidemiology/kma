from ruamel.yaml import YAML

# Create an instance of YAML with the required options
yaml = YAML(typ='unsafe', pure=True)

# Create an ordered dictionary for each section
package = yaml.comments.CommentedMap()
package['name'] = 'kma'

with open('version.h', 'r') as f:
    for line in f:
        if line.startswith('#define'):
            package['version'] = line.split()[2].replace("\"", "")

source = yaml.comments.CommentedMap()
source['url'] = 'https://bitbucket.org/genomicepidemiology/{}/get/{}.tar.gz'.format(package['name'], package['version'])

build = yaml.comments.CommentedMap()
build['number'] = 0
#build['noarch'] = 'generic'

requirements = yaml.comments.CommentedMap()
requirements['build'] = ['make', '{{ compiler(\'c\') }}']
requirements['host'] = ['zlib']
requirements['run'] = ['zlib']

about = yaml.comments.CommentedMap()
about['home'] = 'https://bitbucket.org/genomicepidemiology/kma'
about['summary'] = 'KMA is mapping a method designed to map raw reads directly against redundant databases, in an ultra-fast manner using seed and extend.'
about['license'] = 'Apache-2.0'

extra = yaml.comments.CommentedMap()
identifiers = yaml.comments.CommentedMap()
identifiers['doi'] = '10.1186/s12859-018-2336-6'
extra['identifiers'] = identifiers

# Create a dictionary for the entire YAML content
data = yaml.comments.CommentedMap()
data['package'] = package
data['source'] = source
data['build'] = build
data['requirements'] = requirements
data['about'] = about
data['extra'] = extra

# Serialize the data to a YAML string
yaml_str = yaml.dump(data).replace("\"{{", "{{").replace("}}\"", "}}")

# Print the YAML string
print(yaml_str)

# Write the YAML string to a file
with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
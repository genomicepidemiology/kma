import ruamel.yaml as yaml

# Create an ordered dictionary for each section
package = yaml.comments.CommentedMap()
package['name'] = 'kma'
package['version'] = '1.4.14'

source = yaml.comments.CommentedMap()
source['url'] = 'https://bitbucket.org/genomicepidemiology/kma/get/1.4.14.tar.gz'

build = yaml.comments.CommentedMap()
build['number'] = 0

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

# Serialize the data to YAML and print it
yaml_str = yaml.dump(data, Dumper=yaml.RoundTripDumper).replace("\"{{", "{{").replace("}}\"", "}}")
print(yaml_str)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
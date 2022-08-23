from sys import argv


def parse_names(names_file):
    unique_names = set()
    with open(names_file, 'r') as names:
        for line in names:
            line = line.strip()
            if line:
                substrates = line.split('|')
                for substrate in substrates:
                    unique_names.add(substrate)

    unique_names = list(unique_names)
    unique_names.sort()

    return unique_names


if __name__ == "__main__":
    for name in parse_names(argv[1]):
        print(name)


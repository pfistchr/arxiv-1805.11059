import fractions
import itertools

def fractionfromhex(string):
    if (len(string) != 16) or (string[1] != '.'):
        raise RuntimeError('oops')

    return fractions.Fraction(int(string[0] + string[2:], 16), 1 << 56)

permutations = tuple(itertools.permutations(range(3)))

def computeextremepoints(lower, upper):
    points = set()

    for permutation in permutations:
        lowervalue = lower[permutation[0]]
        uppervalue = upper[permutation[1]]
        othervalue = 1 - lowervalue - uppervalue

        if (othervalue < lower[permutation[2]]) or (othervalue > upper[permutation[2]]):
            raise RuntimeError('oops')

        point = [None]*3
        point[permutation[0]] = lowervalue
        point[permutation[1]] = uppervalue
        point[permutation[2]] = othervalue
        points.add(tuple(point))

    return sorted(points)

infile = open('input.txt', 'r')
outfile = open('input.wl', 'w')
todo = [([fractions.Fraction(0)] * 6, [fractions.Fraction(1)] * 6)]

while todo:
    lower, upper = todo.pop()
    lower[0] = max(lower[0], 1 - upper[1] - upper[2])
    lower[1] = max(lower[1], 1 - upper[0] - upper[2])
    lower[2] = max(lower[2], 1 - upper[0] - upper[1])
    lower[3] = max(lower[3], 1 - upper[4] - upper[5])
    lower[4] = max(lower[4], 1 - upper[3] - upper[5])
    lower[5] = max(lower[5], 1 - upper[3] - upper[4])
    upper[0] = min(upper[0], 1 - lower[1] - lower[2])
    upper[1] = min(upper[1], 1 - lower[0] - lower[2])
    upper[2] = min(upper[2], 1 - lower[0] - lower[1])
    upper[3] = min(upper[3], 1 - lower[4] - lower[5])
    upper[4] = min(upper[4], 1 - lower[3] - lower[5])
    upper[5] = min(upper[5], 1 - lower[3] - lower[4])

    if not all((x >= 0) and (x < y) and (y <= 1) for x, y in zip(lower, upper)):
        raise RuntimeError('oops')

    line = infile.readline()

    if line.startswith('v'):
        tokens = line[1:].split()

        qxextremepoints = computeextremepoints(lower[:3], upper[:3])
        qyextremepoints = computeextremepoints(lower[3:], upper[3:])
        parts = list()

        for qx in qxextremepoints:
            for qy in qyextremepoints:
                qxy = [qx[index // 3] * qy[index % 3] for index in range(9)]
                parts.append(','.join(str(z) for z in qxy))

        outfile.write('verify[{},\n'
                      '{{{{{}}}}},\n'
                      '{{{{{}}}}}];\n'
                      '\n'
                      .format(fractionfromhex(tokens[0]),
                              '},\n{'.join(str(fractionfromhex(z)) for z in tokens[1:10]),
                              '},\n{'.join(parts)))
        continue

    splitindex = {'a\n': 0, 'b\n': 1, 'c\n': 2, 'd\n': 3, 'e\n': 4, 'f\n': 5}[line]
    middlevalue = (lower[splitindex] + upper[splitindex]) / 2

    modifiedlower = lower[:]
    modifiedupper = upper[:]
    modifiedlower[splitindex] = middlevalue
    modifiedupper[splitindex] = middlevalue
    todo.append((modifiedlower, upper))
    todo.append((lower, modifiedupper))

outfile.write('Print["finish"];\n')

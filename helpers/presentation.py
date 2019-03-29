from IPython.core.display import HTML


def show_list(data, mode='bullet'):
    assert mode in {'bullet', 'numeric'}
    output = ''

    if mode == 'bullet':
        container = 'ul'
    else:
        container = 'ol'

    output += f'<{container}>'

    for element in data:
        output += '<li>' + str(element)

    output += f'</{container}>'

    return HTML(output)


def compare_sets(a, b):
    a = set(a)
    b = set(b)

    additional_in_a = a - b
    additional_in_b = b - a

    out = []

    for i, difference in enumerate([additional_in_a, additional_in_b]):
        if difference:
            if len(difference) < 10:
                difference_text = difference
            else:
                dl = list(difference)
                difference_text = (
                    '{'
                        + ', '.join(dl[:3])
                        + ', ..., '
                        + ', '.join(dl[-3:])
                    + '}' + f', {len(difference)} in total'
                )
            out += [f'The {"first" if i == 0 else "second"} set has additional elements: {difference_text}']

    if a == b:
        assert not additional_in_a and not additional_in_b
        out = ['The sets are equal']

    return HTML('<br>'.join(out))
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


def compare_sets(a, b, percentage=None):
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
                dl = sorted(difference)
                total_diff = len(difference)

                if percentage:
                    total = len(a.union(b))
                    total_diff = f'{total_diff} ({total_diff / total * 100:.2f}%)'

                difference_text = (
                    '{'
                        + ', '.join(dl[:3])
                        + ', ..., '
                        + ', '.join(dl[-3:])
                    + '}' + f', {total_diff} in total'
                )

            out += [f'The {"first" if i == 0 else "second"} set has additional elements: {difference_text}']

    if a == b:
        assert not additional_in_a and not additional_in_b
        out = ['The sets are equal']

    return HTML('<br>'.join(out))

#!/usr/bin/env python

import click
import pandas
import plotnine as pn
import math


def set_phase(row):
    phase = row['phase']
    if phase == 'Communication':
        return 'MPI_Irecv'
    elif phase == 'Push_to_buckets':
        return 'Push to Buckets'
    elif phase == 'Waiting_and_Processing':
        return 'Insert into Buckets'
    else:
        return phase

def stackedbars(df):
    df['phase'] = df.apply(lambda row: set_phase(row), axis=1)
    df.phase = pandas.Categorical(df.phase,
                                  categories=['Total', 'MPI_Irecv', 'Push to Buckets', 'Insert into Buckets'], ordered=True)
    p = (pn.ggplot(df, pn.aes(x='m', y='time', fill='phase'))
         + pn.geom_bar(position='stack', stat='identity')
         + pn.scale_x_continuous(trans='log2', labels=lambda ls: ['$2^{' + str(int(math.log2(l))) + '}$' for l in ls])
         + pn.labs(x='m: Number of machines',
                   y='Execution Time (s)',
                   fill='')
         + pn.scale_fill_brewer(type='qual', palette='Set2')
         + pn.theme(legend_position='top',
                    text=pn.element_text(size=12)))
    p += pn.facet_wrap('Component')
    return p

def set_component(row):
    phase = row['phase']
    if phase == 'Total':
        return 'Receiver Total'
    if phase in ['Push_to_buckets', 'Communication']:
        return 'Communicating Thread'
    if phase == 'Waiting_and_Processing':
        return 'Bucketing Thread'


@click.command()
@click.option('--data', help='The path to the file with the data.')
@click.option('--output', help='The output filename for the plot.')
def main(data, output):
    df = pandas.read_csv(data, sep='\t')

    df = df.melt(id_vars=['m'],
                 value_vars=[
                     'Total', 'Waiting_and_Processing', 'Communication',
                     'Push_to_buckets'],
                 var_name='phase', value_name='time')
    df['Component'] = df.apply(lambda row: set_component(row), axis=1)

    df['m'] = df['m'].astype(str).astype(int)

    p = stackedbars(df)
    factor = 1.2
    p.save(output, width=7*factor, height=5*factor)

if __name__ == '__main__':
    main()

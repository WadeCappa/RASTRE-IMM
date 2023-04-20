#!/usr/bin/env python

import click
import pandas
import plotnine as pn
import mizani
import math


def lineplot(df):
    p = (pn.ggplot(df, pn.aes(x='m', y='time', color='phase'))
         + pn.geom_point()
         + pn.geom_line()
         + pn.scale_color_brewer(
             type='qual',
             labels=['Receiver AllToAll', 'Receiver Sampling', 'Receiver Seed Selection',
                     'Sender AllToAll', 'Sender Sampling', 'Sender Seed Selection',
                     'Total'])
         + pn.scale_x_continuous(trans='log2')
         + pn.labs(x='Cluster Nodes (#)',
                   y='Execution Time (s)',
                   color='Phase'))
    return p

def set_phase(row):
    phase = row['phase']
    if phase.startswith('S-') or phase.startswith('R-'):
        return phase[2:]
    else:
        return phase

def stackedbars(df):
    df['phase'] = df.apply(lambda row: set_phase(row), axis=1)
    df.phase = pandas.Categorical(df.phase,
                                   categories=['Total', 'SeedSelect', 'All2All', 'Sampling'], ordered=True)
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

def bars(df):
    data = df.groupby(['m', 'Component']).sum()
    data = data.reset_index()
    p = (pn.ggplot(data, pn.aes(x='m', y='time', fill='Component'))
         + pn.geom_bar(position='dodge', stat='identity')
         + pn.scale_x_continuous(trans='log2')
         + pn.labs(x='Cluster Nodes (#)',
                   y='Execution Time (s)',
                   color=''))
    return p

def set_component(row):
    if row['phase'] == 'Total':
        return 'Total'
    if row['phase'].startswith('R'):
        return 'Receiver'
    if row['phase'].startswith('S'):
        return 'Sender'


@click.command()
@click.option('--data', help='The path to the file with the data.')
@click.option('--output', help='The output filename for the plot.')
def main(data, output):
    df = pandas.read_csv(data, sep='\t')

    df = df.melt(id_vars=['m'],
                 value_vars=[
                     'Total', 'S-Sampling', 'S-All2All', 'S-SeedSelect',
                     'R-Sampling', 'R-All2All', 'R-SeedSelect'],
                 var_name='phase', value_name='time')
    df['Component'] = df.apply(lambda row: set_component(row), axis=1)
    df['m'] = df['m'].astype(str).astype(int)
    p = stackedbars(df)
    factor = 1.2
    p.save(output, width=7*factor, height=5*factor)

if __name__ == '__main__':
    main()

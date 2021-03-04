#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2020
# File   : signature_frame.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 04.06.2020

from pyfiglet import Figlet

from goetia.cli.cli import format_filenames, TextBlock


class SignatureStreamFrame:

    def __init__(self, term, args, component_name='sourmash stream'):
        try:
            import plotille as pt
        except ImportError:
            print('plotille is required for terminal plotting: '
                  'install with `pip install plotille`',
                  file=sys.stderr)
            sys.exit(1)

        self.term = term
        self.args   = args

        self.distance_figure = pt.Figure()
        self.distance_figure.width = 60
        self.distance_figure.height = 20
        self.distance_figure.x_label = '# sequences'
        self.distance_figure.y_label = 'Δdistance'
        self.distance_figure.set_x_limits(min_=0)
        self.distance_figure.set_y_limits(min_=0)
        self.distance_figure.color_mode = 'byte'
        self.distance_figure.origin = False

        self.hist_figure = pt.Figure()
        self.hist_figure.width = 30
        self.hist_figure.height = 20
        self.hist_figure.x_label = 'Δdistance'
        self.hist_figure.y_label = 'counts'
        self.hist_figure.color_mode = 'byte'

        name_block = Figlet(font='doom').renderText('goetia').rstrip().split('\n')
        name_block[-1] = name_block[-1] + f' {__version__}'
        name_block.append(f'\n{term.bold}{component_name}{term.normal}')
        self.name_block = TextBlock('\n'.join(name_block))

        param_block = term.normal + term.underline + (' ' * 40) + term.normal + '\n\n'
        param_block += f'{term.bold}window size:       {term.normal}{args.window_size}\n'\
                       f'{term.bold}tick length:       {term.normal}{args.tick_length}\n'\
                       f'{term.bold}distance metric:   {term.normal}{args.distance_metric}\n'\
                       f'{term.bold}saturation policy: {term.normal}{args.saturation_policy}\n'\
                       f'{term.bold}cutoff:            {term.normal}{args.cutoff}'

        self.param_block = TextBlock(param_block)
        metric_text = term.normal + term.underline + (' ' * 40) + term.normal + '\n'
        metric_text += '{term.move_down}{term.bold}'
        metric_text += '# sequences:  '.ljust(20)
        metric_text += '{term.normal}{n_reads:,}\n'
        metric_text += '{term.bold}'
        metric_text += 'Δdistance:    '.ljust(20)
        metric_text += '{term.normal}{prev_d:.6f}\n'
        metric_text += '{term.bold}'
        metric_text += 'windowed(Δdistance): '.ljust(20)
        metric_text += '{term.normal}{stat:.6f}\n'
        metric_text += '{term.underline}' + (' ' * 40) + '{term.normal}'
        self.metric_text = metric_text

        #message_text = term.normal + term.underline + (' ' * 40) + term.normal + '\n'
        message_text = '{messages}'
        self.message_text = message_text

        self.distances = []
        self.distances_t = []

    def __del__(self):
        self._print(self.term.move_down * (self.distance_figure.height + 10))

    def _print(self, *args, **kwargs):
        self.term.stream.write(*args, **kwargs)

    def metric_block(self, n_reads, prev_d, stat):
        return TextBlock(self.metric_text.format(n_reads=n_reads, prev_d=prev_d, stat=stat, term=self.term))

    def message_block(self, messages):
        text = self.message_text.format(messages='\n'.join(messages))
        return TextBlock(text)

    def figure_block(self, distances, distances_t):
        self.distance_figure.clear()
        self.distance_figure.plot(distances_t, distances, lc=10, interp=None)
        text = self.distance_figure.show()
        return TextBlock(text)

    def hist_block(self, distances):
        self.hist_figure.clear()
        distances = np.log(distances)
        self.hist_figure.histogram(distances, lc=20)
        text = self.hist_figure.show()
        return TextBlock(text)

    async def draw(self, n_reads=0, prev_d=1.0, stat=1.0,
                   messages=None, draw_dist_plot=True,
                   draw_dist_hist=False):

        term = self.term
        self.distances.append(prev_d)
        self.distances_t.append(n_reads)

        figure_thr = await curio.spawn_thread(self.figure_block, self.distances, self.distances_t)
        metrics_thr = await curio.spawn_thread(self.metric_block, n_reads, prev_d, stat)
        
        figure = await figure_thr.join()
        metrics = await metrics_thr.join()

        term.stream.write(term.clear_eos)

        if draw_dist_plot:
            with term.location():
                figure.draw(term, max(self.name_block.width, self.param_block.width) + 1, 1)

        with term.location():
            self.name_block.draw(term, 0, shift_y=1)

        with term.location():
            self.param_block.draw(term, 0, shift_y=self.name_block.height+4)

        with term.location():
            metrics.draw(term,
                         0,
                         self.name_block.height + self.param_block.height + 4)

        if messages is not None:
            with term.location():
                ypos = _metric_block.height + self.name_block.height + self.param_block.height + 9
                self.message_block(messages).draw(term, 0, ypos)

        if draw_dist_hist and distances is not None:
            with term.location():
                xpos = max(self.name_block.width, self.param_block.width) + 2
                if draw_dist_plot:
                    xpos += self.distance_figure.width + 5
                self.hist_block(distances).draw(term, xpos, 1)


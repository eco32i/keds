from ggplot import *

def plot_shape(temp):
    df = df_dHSR1[(df_dHSR1['Temp'] == float(temp)) & (df_dHSR1['pos'] < 639)]
    p = ggplot(df, aes(xmin='pos', xmax='pos+1', ymin=0, ymax='theta28')) \
        + geom_rect() \
        + scale_x_continuous(name='nt position', limits=(-10,610)) \
        + scale_y_continuous(name='SHAPE reactivity')
    print p


def plot_diff(temp1, temp2):
    colt1 = 'diff_{:.1f}'.format(float(temp1))
    colt2 = 'diff_{:.1f}'.format(float(temp2))
    df['diff'] = df[colt2] - df[colt1]
    p = ggplot(df, aes(xmin='pos', xmax='pos+1', ymin=0, ymax='diff')) \
        + geom_rect() \
        + scale_x_continuous(name='nt position', limits=(-10,650), breaks=(100,200,300,400,500,600)) \
        + scale_y_continuous(name='SHAPE reactivity diff') \
        + theme(axis_text_x=element_text(size=20),
                axis_title_x=element_text(size=20),
                axis_text_y=element_text(size=20),
                axis_title_y=element_text(size=20))
    print p


def SHAPE_map(df, temp, field='theta28', num_cols=60):
    data = shape_df(df, temp, field=field)
    num_rows = data.shape[0] / num_cols
    if data.shape[0] % num_cols:
        num_rows += 1

    data['row'] = [(i / num_cols) + 1 for i in range(data.shape[0])]
    data['col'] = [(i - i / num_cols * num_cols) + 1 for i in range(data.shape[0])]

    row_range = [str(x) for x in reversed(range(1, num_rows + 1))]
    col_range = [str(x) for x in range(1, num_cols + 1)]
    max_theta = data[field].max()

    alpha = data['theta28'] / max_theta

    source = ColumnDataSource(
        data=dict(
            row=[str(x) for x in data['row']],
            col=[str(x) for x in data['col']],
            base=data['base'],
            theta=data[field],
            alphas=alpha,
        )
    )

    p = figure(title='Drosophila HSR1 SHAPE reactivity at {:.1f} deg'.format(temp),
               background_fill='#fdf6e3',
               tools='resize,hover,save',
               x_range=col_range, y_range=row_range, plot_width=1000)

    p.rect('col', 'row', 0.9, 0.9, source=source, alpha='alphas', color='#dc322f', line_color=None)

    text_props = {
        'source': source,
        'angle': 0,
        'color': '#586e75',
        'text_align': 'center',
        'text_baseline': 'middle',
    }

    p.text(x=dict(field='col', units='data'),
           y=dict(field='row', units='data'),
           text=dict(field='base', units='data'),
           text_font_size='10pt', text_font_style='bold', **text_props)

    p.grid.grid_line_color = None
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [
        ('base', '@base'),
        ('theta', '@theta')
    ]
    show(p)


def SHAPE_diffmap(df, temp1, temp2, notitle=False, field='theta28', num_cols=60):
    from bokeh.models import LinearAxis, SingleIntervalTicker
    data = diff_df(df_plus, temp1, temp2, field=field)
    num_rows = data.shape[0] / num_cols
    if data.shape[0] % num_cols:
        num_rows += 1

    data['row'] = [(i / num_cols) + 1 for i in range(data.shape[0])]
    data['col'] = [(i - i / num_cols * num_cols) + 1 for i in range(data.shape[0])]
    data['color'] = '#dc322f'
    data.loc[data['diff'] < 0, 'color'] = '#268bd2'

    row_range = [str(x) for x in reversed(range(1, num_rows + 1))]
    col_range = [str(x) for x in range(1, num_cols + 1)]

    max_diff = max(data['diff'].max(), abs(data['diff'].min())) or 1

    data['alpha'] = abs(data['diff'] / max_diff)

    source = ColumnDataSource(
        data=dict(
            row=[str(x) for x in data['row']],
            col=[str(x) for x in data['col']],
            base=data['base'],
            color=data['color'],
            alphas=data['alpha'],
        )
    )

    if notitle:
        title = None
    else:
        title='Drosophila HSR1 SHAPE reactivity difference {:.1f} vs {:.1f} deg'.format(temp2, temp1)
    p = figure(title=title,
               background_fill='#fdf6e3',
               tools='resize,hover,save', x_axis_type=None,
               x_range=col_range, y_range=row_range, plot_width=1000)

    p.rect('col', 'row', 0.9, 0.9, source=source, alpha='alphas', color='color', line_color=None)

    text_props = {
        'source': source,
        'angle': 0,
        'color': '#586e75',
        'text_align': 'center',
        'text_baseline': 'middle',
    }

    p.text(x=dict(field='col', units='data'),
           y=dict(field='row', units='data'),
           text=dict(field='base', units='data'),
           text_font_size='10pt', text_font_style='bold', **text_props)

    p.grid.grid_line_color = None
    ticker = SingleIntervalTicker(interval=10)
    xaxis = LinearAxis(ticker=ticker)
    p.add_layout(xaxis, 'below')

    p.axis.major_label_text_font_size = "22pt"

    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [
        ('base', '@base'),
        ('diff', '@alphas')
    ]
    show(p)


def SHAPE_range(df, temps, field='theta28', title=None, norm_local=True):
    from bokeh.models import LinearAxis, SingleIntervalTicker
    data = df[(df['pos'] > 0) & (df['pos'] < 605)][['base', 'pos', 'Temp', field]].copy()
    row_range = [str(float(x)) for x in sorted(temps)]
    col_range = [str(x) for x in set(data['pos'])]
    if norm_local:
        data['alphas'] = 0
        # We want to compute alphas values separately for each experiment (Temp)
        for temp in temps:
            max_theta = data[data['Temp'] == temp][field].max()
            data.loc[data['Temp'] == temp, ['alphas']] = data[data['Temp'] == temp][field] / max_theta
    else:
        max_theta = data[field].max()
        data['alphas'] = data[field] / max_theta

    source = ColumnDataSource(
        data=dict(
            temp=[str(float(x)) for x in data['Temp']],
            col=[str(x) for x in data['pos']],
            base=data['base'],
            theta=data[field],
            alphas=data['alphas'],
        )
    )

    p = figure(title=title, background_fill='#fdf6e3',
               tools='resize,save,reset,box_zoom', x_axis_type=None,
               x_range=col_range, y_range=row_range, plot_width=1000)

    p.rect('col', 'temp', 0.9, 0.9, source=source, alpha='alphas', color='#dc322f', line_color=None)
    p.grid.grid_line_color = None
    ticker = SingleIntervalTicker(interval=100, num_minor_ticks=10)
    xaxis = LinearAxis(ticker=ticker)
    p.add_layout(xaxis, 'below')
    show(p)


def SHAPE_diffrange(df, temps, basetemp=36.9, title=None, field='theta28', norm_local=True):
    from bokeh.models import LinearAxis, SingleIntervalTicker
    data = df[(df['pos'] > 0) & (df['pos'] < 605)][['base', 'pos', 'Temp', field]].copy()
    row_range = [str(float(x)) for x in sorted(temps)]
    col_range = [str(x) for x in set(data['pos'])]

    data['diff'] = 0
    for temp in temps:
        data.loc[data['Temp'] == temp, ['diff']] = data[data['Temp'] == temp][field] - data[data['Temp'] == basetemp][field]

    data['color'] = '#dc322f'
    data.loc[data['diff'] < 0, 'color'] = '#268bd2'

    if norm_local:
        data['alphas'] = 0
        # We want to compute alphas values separately for each experiment (Temp)
        for temp in temps:
            max_diff = max(data[data['Temp'] == temp]['diff'].max(), abs(data[data['Temp'] == temp]['diff'].min())) or 1
            data.loc[data['Temp'] == temp, ['alphas']] = data[data['Temp'] == temp]['diff'] / max_diff
    else:
        max_diff = max(data['diff'].max(), abs(data['diff'].min())) or 1
        data['alphas'] = data['diff'] / max_diff

    data.loc[data['Temp'] == basetemp, 'alphas'] = 0.15
    source = ColumnDataSource(
        data=dict(
            temp=[str(float(x)) for x in data['Temp']],
            col=[str(x) for x in data['pos']],
            base=data['base'],
            theta=data[field],
            alphas=data['alphas'],
            color=data['color']
        )
    )

    p = figure(title=title,  background_fill='#fdf6e3',
               tools='resize,save,reset,box_zoom', x_axis_type=None,
               x_range=col_range, y_range=row_range, plot_width=1000)

    p.rect('col', 'temp', 0.9, 0.9, source=source, alpha='alphas', color='color', line_color=None)
    p.grid.grid_line_color = None
    ticker = SingleIntervalTicker(interval=100, num_minor_ticks=10)
    xaxis = LinearAxis(ticker=ticker)
    p.add_layout(xaxis, 'below')
    p.axis.major_label_text_font_size = "22pt"

    show(p)


def SHAPE_change(df, temps, field='theta28', threshold=3, local=False):
    data = df[(df['pos'] > 0) & (df['pos'] < 605)][['base', 'pos', 'Temp', field]].copy()

    data['diff'] = 0
    basetemp = min(temps)

    for temp in temps:
        data.loc[data['Temp'] == temp, ['diff']] = abs(data[data['Temp'] == temp][field] - data[data['Temp'] == basetemp][field])

    if not local:
        mean_diff = data['diff'].mean()
        sigma = data['diff'].std()
        #print "Mean:\t%f\tsigma:\t%f" % (mean_diff,sigma)

    _df = pd.DataFrame()
    _df['Temp'] = temps[1:]
    _df['sum'] = 0
    for temp in temps:
        if temp != basetemp:
            if local:
                mean_diff = data[data['Temp'] == temp]['diff'].mean()
                sigma = data[data['Temp'] == temp]['diff'].std()
            tdiff = data[(data['Temp'] == temp) & (data['diff'] - mean_diff > threshold * sigma)]['diff']
            _df.loc[_df['Temp'] == temp, 'sum'] = tdiff.sum() / tdiff.shape[0]
    return _df

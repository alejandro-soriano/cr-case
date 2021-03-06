import os

import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import itertools

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq

import plotly.graph_objs as go

import sections as sc


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

server = app.server


#Define dummy input parameters
V_y = 100 #N
V_x = 100. #N

t_w = 7.1 #mm
t_f_u = 10.7 #mm
t_f_l = 10.7 #mm
h = 300. #mm
w_u = 300. #mm
w_l = 300. #mm
num = 7 #Initial dumber of discretization elements
disc = 'N'

data, layout = sc.gen_data_i_beam(V_y,V_x,t_w,t_f_u,t_f_l,h,w_u,w_l,num,disc)

# Generate a Plotly figure object that is used in conjunction with Dash widgets
fig = go.Figure(data= data, layout = layout)

#Options provided in dropdown
opts = [{'label': 'I Shape', 'value': 'I'},
        {'label': 'T Shape', 'value': 'T'},
        {'label': 'C Shape', 'value': 'C'},
        {'label': 'U Shape', 'value': 'U'},
        {'label': 'Rect. Shape', 'value': 'S'},
        # {'label': 'Y Shape', 'value': 'Y'},
        {'label': 'Circ. w/ slit', 'value': 'O'},]

opts_disc = [{'label': 'Yes', 'value': 'Y'},
        {'label': 'No', 'value': 'N'},]

#Definition of the layout of the application. Includes the main figure and the widgets that customize the input
drop_section = html.Div([
        html.Label("Choose cross-section type"),
        #Dropdown definition
        dcc.Dropdown(id = 'opt', options = opts,
                    value = 'I')
            ], style = {'width': '168px',
                        'fontSize' : '12px',
                        'padding-left' : '100px',
                        'display': 'table-cell'})

drop_disc = html.Div([
        html.Label("Show discretization?"),
        #Dropdown definition
        dcc.Dropdown(id = 'opt_disc', options = opts_disc,
                    value = 'N')
            ], style = {'width': '168px',
                        'fontSize' : '12px',
                        'padding-left' : '20px',
                        'display': 'table-cell'})

input_V =    html.Div([ html.Label("Define shear load (V) [N]"),
    #Input textbox
    dcc.Input(
        id='inp_V',
        placeholder='Enter a value...',
        type='number',
        value=V_y)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_h = html.Div([
    html.Label("Define section height (h) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_h',
        placeholder='Enter a value...',
        type='number',
        value=h)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '100px',
                    'display': 'table-cell'})

input_h_l = html.Div([
    html.Label("Define left height (h1) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_h_l',
        placeholder='Enter a value...',
        type='number',
        value=h)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '100px',
                    'display': 'table-cell'})

input_h_r = html.Div([
    html.Label("Define right height (h2) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_h_r',
        placeholder='Enter a value...',
        type='number',
        value=h/2.)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_w1 = html.Div([
    html.Label("Define top width (w1) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_w1',
        placeholder='Enter a value...',
        type='number',
        value=w_u)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})
input_w = html.Div([
    html.Label("Define section width (w) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_w',
        placeholder='Enter a value...',
        type='number',
        value=w_u)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_w2 = html.Div([
    html.Label("Define bottom width (w2) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_w2',
        placeholder='Enter a value...',
        type='number',
        value=w_l)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_t1 =  html.Div([
    html.Label("Define top thickness (t1) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_t1',
        placeholder='Enter a value...',
        type='number',
        value=t_f_u)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_t2 = html.Div([
    html.Label("Define bottom thick. (t2) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_t2',
        placeholder='Enter a value...',
        type='number',
        value=t_f_l)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_t3 = html.Div([
    html.Label("Define web thickness (t3) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_t3',
        placeholder='Enter a value...',
        type='number',
        value=t_w)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})
input_t3_l = html.Div([
    html.Label("Define left web thick. (t3) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_t3_l',
        placeholder='Enter a value...',
        type='number',
        value=t_w)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_t4_r = html.Div([
    html.Label("Define right w. thick. (t4) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_t4_r',
        placeholder='Enter a value...',
        type='number',
        value=t_w)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_t5 =  html.Div([
    html.Label("Define circle thick. (t5) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_t5',
        placeholder='Enter a value...',
        type='number',
        value=t_f_u)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'})

input_rad =  html.Div([
    html.Label("Define circle radius (rad) [mm]"),
    #Input textbox
    dcc.Input(
        id='inp_rad',
        placeholder='Enter a value...',
        type='number',
        value=w_u/2.)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '100px',
                    'display': 'table-cell'})


#Figure
graph_fig = dcc.Graph(
        id='flow',
        figure=fig
    )

input_slider = html.Div([
    html.Label('Define discretization'),
    #Slider object
    dcc.Slider(id='slider',min=3,max=21,step=2,value=7,
    marks={int(i): '{}'.format(i) if i == 3 else str(i) for i in np.arange(3, 23,2)},
    ),],
        style = {'width' : '80%',
                'fontSize' : '12px',
                'padding-left' : '100px',
                'display': 'inline-block'})

output_slider = html.Div(id='slider-output-container')


app.layout = html.Div([
    html.Div([html.H4("Shear Flow App")], style={'font-size': '12px', 'textAlign': "left", "padding-left": '100px'}),
    html.Div(["Calculates and draws the shear flow field for various engineering cross-sections. Having chosen a section \
        the shear flow is calculated according to its linear exact solution in the continous form. An approximation based \
        on a defined discretization level is presented to the side. The app supports symmetric and non-symmetric cross-sections. \
        More information ",html.A('here', href='https://github.com/alejandro-soriano/cr-case'), "."], style={'width': '85%', 'textAlign': "left", "padding-left": '100px'}),
    html.Div([html.Hr()], style={"margin-before": '10px'}),
    html.Div([drop_section,drop_disc,input_V,]),
    html.Div(
        id='controls-container'
    ),
    html.Div(
        id='output-container'
    ),
    html.Div(
        id='slider-container'
    )
    # input_slider,
    # output_slider
])

def generate_control_id(value):
    return 'Control {}'.format(value)

DYNAMIC_CONTROLS = {
    None: html.Div([input_h, input_w1, input_w2, input_t1, input_t2,input_t3
        ]), 
    'I': html.Div([input_h, input_w1, input_w2, input_t1, input_t2,input_t3
        ]),
    'T': html.Div([input_h, input_w1, input_t1, input_t3
        ]),
    'C': html.Div([input_h, input_w1, input_w2, input_t1, input_t2,input_t3
        ]),
    'S': html.Div([input_h, input_w, input_t1, input_t2, input_t3_l, input_t4_r
        ]),
    'U': html.Div([input_h_l, input_h_r, input_w, input_t2, input_t3_l, input_t4_r
        ]),
    'O': html.Div([input_rad, input_t5
        ])
}

DYNAMIC_CONTROLS_SLIDER = {
    None: html.Div([
        ]), 
    'Y': html.Div([input_slider,
                   output_slider
        ]),
    'N': html.Div([
        ])
}

@app.callback(
    dash.dependencies.Output('controls-container', 'children'),
    [dash.dependencies.Input('opt', 'value'),
    dash.dependencies.Input('opt_disc', 'value')])
def display_controls(datasource_1_value,datasource_2_value):
    # generate dynamic controls based off of the cross section source selections
    return html.Div([
        DYNAMIC_CONTROLS[datasource_1_value]
    ])

@app.callback(
    dash.dependencies.Output('slider-container', 'children'),
    [dash.dependencies.Input('opt', 'value'),
    dash.dependencies.Input('opt_disc', 'value')])
def display_controls(datasource_1_value,datasource_2_value):
    # generate dynamic controls based off of the cross section source selections
    return html.Div([
        DYNAMIC_CONTROLS_SLIDER[datasource_2_value]
    ])

def generate_output_id(value1,value2):
    return '{} {} container'.format(value1,value2)

@app.callback(
    dash.dependencies.Output('output-container', 'children'),
    [dash.dependencies.Input('opt', 'value'),
     dash.dependencies.Input('opt_disc', 'value')])

def display_controls(datasource_1_value,datasource_2_value):
    # create a unique output container for each pair of dyanmic controls
    return html.Div(id=generate_output_id(
        datasource_1_value,
        datasource_2_value
    ))

def generate_output_callback(datasource_1_value,datasource_2_value):
    if datasource_1_value is None:
        vacio = []
        #do nothing
    elif datasource_1_value == 'I':
        if datasource_2_value == 'Y':
            def output_callback(input1,input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_i_beam(input3,V_x,input9,input7,input8,input4,input5,input6,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )
                return html.Div(graph_fig)

        elif datasource_2_value == 'N':
            input1 = 7
            def output_callback(input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_i_beam(input3,V_x,input9,input7,input8,input4,input5,input6,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )
                return html.Div(graph_fig)


    elif datasource_1_value == 'T':
        if datasource_2_value == 'Y':
            def output_callback(input1,input2,input3,input4,input5,input7,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_t_beam(input3,V_x,input9,input7,input4,input5,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
        elif datasource_2_value == 'N':
            input1 = 7
            def output_callback(input2,input3,input4,input5,input7,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_t_beam(input3,V_x,input9,input7,input4,input5,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)


    elif datasource_1_value == 'C':
        if datasource_2_value == 'Y':
            def output_callback(input1,input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_c_beam(input3,V_x,input9,input7,input8,input4,input5,input6,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
        elif datasource_2_value == 'N':
            input1 = 7
            def output_callback(input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_c_beam(input3,V_x,input9,input7,input8,input4,input5,input6,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
    
    elif datasource_1_value == 'S':
        if datasource_2_value == 'Y':
            def output_callback(input1,input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_s_beam(input3,V_x,input8,input9,input6,input7,input4,input5,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
        elif datasource_2_value == 'N':
            input1 = 7
            def output_callback(input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_s_beam(input3,V_x,input8,input9,input6,input7,input4,input5,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
    
    elif datasource_1_value == 'U':
        if datasource_2_value == 'Y':
            def output_callback(input1,input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_u_beam(input3,V_x,input7,input8,input9,input4,input5,input6,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
        elif datasource_2_value == 'N':
            input1 = 7
            def output_callback(input2,input3,input4,input5,input6,input7,input8,input9,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_u_beam(input3,V_x,input7,input8,input9,input4,input5,input6,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)

    elif datasource_1_value == 'O':
        if datasource_2_value == 'Y':
            def output_callback(input1,input2,input3,input4,input5,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_o_beam(input3,V_x,input4,input5,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)
        elif datasource_2_value == 'N':
            input1 = 7
            def output_callback(input2,input3,input4,input5,input10):
                # This function can display different outputs depending on
                # the values of the dynamic controls
                data, layout = sc.gen_data_o_beam(input3,V_x,input4,input5,input1,input10)

                fig = go.Figure(data = data, layout = layout)

                graph_fig = dcc.Graph(
                            id='flow',
                            figure=fig
                        )

                return html.Div(graph_fig)

    return output_callback

app.config.supress_callback_exceptions = True

# create a callback for all possible combinations of dynamic controls
# each unique dynamic control pairing is linked to a dynamic output component
for value1,value2 in itertools.product(
    [o['value'] for o in app.layout['opt'].options],
    [o['value'] for o in app.layout['opt_disc'].options]):
    if value1 == 'I':
        if value2 == 'Y':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w1', 'value'),
                dash.dependencies.Input('inp_w2', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
        elif value2 == 'N':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [#dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w1', 'value'),
                dash.dependencies.Input('inp_w2', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
    elif value1 == 'T':
        if value2 == 'Y':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w1', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t3', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
        elif value2 == 'N':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [#dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w1', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t3', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
    elif value1 == 'C':
        if value2 == 'Y':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w1', 'value'),
                dash.dependencies.Input('inp_w2', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
        elif value2 == 'N':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [#dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w1', 'value'),
                dash.dependencies.Input('inp_w2', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
    elif value1 == 'S':
        if value2 == 'Y':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3_l', 'value'),
                dash.dependencies.Input('inp_t4_r', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
        elif value2 == 'N':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [#dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h', 'value'),
                dash.dependencies.Input('inp_w', 'value'),
                dash.dependencies.Input('inp_t1', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3_l', 'value'),
                dash.dependencies.Input('inp_t4_r', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )    
    elif value1 == 'U':
        if value2 == 'Y':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h_l', 'value'),
                dash.dependencies.Input('inp_h_r', 'value'),
                dash.dependencies.Input('inp_w', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3_l', 'value'),
                dash.dependencies.Input('inp_t4_r', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
        elif value2 == 'N':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [#dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_h_l', 'value'),
                dash.dependencies.Input('inp_h_r', 'value'),
                dash.dependencies.Input('inp_w', 'value'),
                dash.dependencies.Input('inp_t2', 'value'),
                dash.dependencies.Input('inp_t3_l', 'value'),
                dash.dependencies.Input('inp_t4_r', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
    elif value1 == 'O':
        if value2 == 'Y':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_rad', 'value'),
                dash.dependencies.Input('inp_t5', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )
        elif value2 == 'N':
            app.callback(
                dash.dependencies.Output(generate_output_id(value1,value2), 'children'),
                [#dash.dependencies.Input('slider', 'value'),
                dash.dependencies.Input('opt', 'value'),
                dash.dependencies.Input('inp_V', 'value'),
                dash.dependencies.Input('inp_rad', 'value'),
                dash.dependencies.Input('inp_t5', 'value'),
                dash.dependencies.Input('opt_disc', 'value')
                ])(
                generate_output_callback(value1,value2)
            )




#Run application
if __name__ == '__main__':
    app.run_server(debug=True)
import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate

import itertools
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import pandas as pd
import plotly.graph_objs as go

#Define equation of a parabola for shear flow field plotting
def eq_parab(p,e1,e2,vert): #Function used to calculate the parabola constants
    a, b, c = p
    eq1 = -e1[0] + a*e1[1]**2 + b*e1[1] + c
    eq2 = -e2[0] + a*e2[1]**2 + b*e2[1] + c
    eq3 = -vert[0] + a*vert[1]**2 + b*vert[1] + c
    return (eq1,eq2,eq3)


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

df = pd.read_csv(
    'https://gist.githubusercontent.com/chriddyp/' +
    '5d1ea79569ed194d432e56108a04d188/raw/' +
    'a9f9e8076b837d541398e999dcbac2b2826a81f8/'+
    'gdp-life-exp-2007.csv')

#Define input parameters
V_y = 100 #N
V_x = 0. #N

t_w = 6. #mm
t_f_u = 8. #mm
t_f_l = 8. #mm
h = 100. #mm
w_u = 160. #mm
w_l = 50. #mm

num = 15
def gen_data_i_beam(V_y,V_x,t_w,t_f_u,t_f_l,h,w_u,w_l,num):

    #Calculate geometrical properties of section
    A_tot = w_u*t_f_u + t_w*(h - t_f_u - t_f_l) + w_l*t_f_l #mm2
    y_cdg = (w_u*t_f_u*(h - t_f_u/2.) + t_w*(h - t_f_u - t_f_l)*h/2. + w_l*t_f_l**2/2)/A_tot #mm w.r.t. bottom line
    I_x = 1./12.*w_u*t_f_u**3 + 1./12.*w_l*t_f_l**3 + w_u*t_f_u*(h/2. - t_f_u/2.)**2 + w_l*t_f_l*(h/2. - t_f_l/2.)**2 + 1./12.*t_w*(h - t_f_u - t_f_l)**3 #mm4
    I_y = 1./12.*t_f_u*w_u**3 + 1./12.*t_f_l*w_l**3 + 1./12.*(h-t_f_u-t_f_l)*t_w**3
    I_xy = 0.


    #Compute maximum shear flows on section
    q_flange_u = V_y*(t_f_u*w_u/2.*(h/2. - t_f_u/2.))/I_x #N/mm
    q_flange_l = V_y*(t_f_l*w_l/2.*(h/2. - t_f_l/2.))/I_x #N/mm

    Q_w_u = t_f_u*w_u*(h - y_cdg - t_f_u/2.) + (h - y_cdg - t_f_u)**2/2*t_w #mm3
    Q_w_l = t_f_l*w_l*(y_cdg - t_f_l/2.) + (y_cdg - t_f_l)**2/2*t_w #mm3

    q_web = V_y*Q_w_u/I_x #N/mm @ cog




    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    origin_x = - 1.25*max([w_u,w_l,h/2.])/2.
    origin_y = 0.

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1 = go.Scatter(
        x=[origin_x - w_l/2., origin_x + w_l/2. ,origin_x+w_l/2., origin_x + t_w/2.,origin_x + t_w/2.,origin_x+w_u/2.,origin_x+w_u/2.,origin_x - w_u/2.,origin_x - w_u/2., origin_x-t_w/2.,origin_x-t_w/2.,origin_x - w_l/2.,origin_x - w_l/2.],
        y=[origin_y,origin_y,origin_y+t_f_l,origin_y+t_f_l,origin_y+h-t_f_u,origin_y+h-t_f_u,origin_y+h,origin_y+h,origin_y+h-t_f_u,origin_y+h-t_f_u,origin_y+t_f_l,origin_y+t_f_l,origin_y],
        mode='lines',
        name=r'I Beam <br> V = {0:.2f} N <br> h = {1:.0f} mm <br> w1 = {2:.0f} mm <br> w2 = {3:.0f} mm <br> t1 = {4:.0f} mm <br> t2 = {5:.0f} mm <br> t3 = {6:.0f} mm <br> A = {7:.0f} mm2 <br> cog = {8:.2f} mm <br> I = {9:.0f} mm4'.format(V_y,h,w_u,w_l,t_f_u,t_f_l,t_w,A_tot,y_cdg,I_x),
        hoverinfo='text',
        hoveron='fills',
        fill='toself', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    #Entity that plots the location of the shear load application
    trace1_sc = go.Scatter(
        x = [origin_x],
        y = [y_cdg],
        name = 'Shear Center',
        hovertext = 'Shear Center',
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'black', size = 2*t_w),
        visible=[True if c is 0 else False][0])

    #Determination of the shear flow fields in the upper and lower flanges of the I-section
    y_flange_u = np.linspace(0.,q_flange_u,11).tolist() + np.linspace(q_flange_u,0.,11).tolist()
    y_flange_l = np.linspace(0.,q_flange_l,11).tolist() + np.linspace(q_flange_l,0.,11).tolist()

    #Entity that plots the upper flange shear flow
    trace2 = go.Scatter(
        x= np.linspace(origin_x - w_u/2.,origin_x,11).tolist() + np.linspace(origin_x,origin_x + w_u/2.,11).tolist(),
        y= np.linspace(origin_y + h,origin_y + h + h/14.*q_flange_u/max([q_flange_u,q_flange_l]),11).tolist() + np.linspace(origin_y + h + h/14.*q_flange_u/max([q_flange_u,q_flange_l]),origin_y + h,11).tolist(),
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in y_flange_u],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Entity that plots the lower flange shear flow
    trace3 = go.Scatter(
        x= np.linspace(origin_x - w_l/2.,origin_x,11).tolist() + np.linspace(origin_x,origin_x + w_l/2.,11).tolist(),
        y= np.linspace(origin_y ,origin_y - h/14.*q_flange_l/max([q_flange_u,q_flange_l]),11).tolist() + np.linspace(origin_y - h/14.*q_flange_l/max([q_flange_u,q_flange_l]),origin_y ,11).tolist(),
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in y_flange_l],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Computation of parabola constats for plotting the shear flow on web
    #Points through which the parabola-shaped shear flow in the web goes
    e1_coord = [origin_x - max([w_u/2.,w_l/2.]) - h/14.*q_flange_l/max([q_flange_u,q_flange_l]),t_f_l/2.]
    e2_coord = [origin_x - max([w_u/2.,w_l/2.])- h/14.*q_flange_u/max([q_flange_u,q_flange_l]),h - t_f_u/2.]
    vert_coord = [origin_x -max([w_u/2.,w_l/2.])- h/14.*q_web/max([q_flange_u,q_flange_l]),y_cdg]

    #Solve for the parabola constants
    a_coord, b_coord, c_coord = opt.fsolve(eq_parab,[1/180.,-7./18.,713./36.],args=(e1_coord,e2_coord,vert_coord))

    y_par = np.linspace(origin_y+t_f_l/2.,origin_y+h-t_f_u/2.,21)
    x_par_coord = a_coord*(y_par)**2 + b_coord*y_par + c_coord

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par = y_par.tolist()
    y_par.append(origin_y+h-t_f_u/2.)
    y_par.insert(0,origin_y+t_f_l/2.)

    x_par_coord = x_par_coord.tolist()
    x_par_coord.append(origin_x - max([w_u/2.,w_l/2.]))
    x_par_coord.insert(0,origin_x - max([w_u/2.,w_l/2.]))

    #Points through which the parabola-shaped shear flow in the web goes in order to define the appropriate tags
    e1_flow = [q_flange_l,t_f_l/2.]
    e2_flow = [q_flange_u,h - t_f_u/2.]
    vert_flow = [q_web, y_cdg]

    #Solve for the parabola constants
    a_flow, b_flow, c_flow = opt.fsolve(eq_parab,[1/180.,-7./18.,713./36.],args=(e1_flow,e2_flow,vert_flow))

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow = np.linspace(origin_y+t_f_l/2.,origin_y+h-t_f_u/2.,21)
    x_par_flow = a_flow*(y_par_flow)**2 + b_flow*y_par_flow + c_flow

    x_par_flow = x_par_flow.tolist()
    x_par_flow.append(q_flange_u)
    x_par_flow.insert(0,q_flange_l)

    #Entity that plots the web shear flow
    trace4 = go.Scatter(
        x=x_par_coord,
        y=y_par,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in x_par_flow],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Set the origin of coordinates for the right-most figure (idealized)
    origin_xd = 1.25*max([w_u,w_l,3*h/4.])/2.
    origin_yd = 0.

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = np.linspace(origin_xd - w_u/2. ,origin_xd+w_u/2.,num).tolist()
    x_fl = np.linspace(origin_xd - w_l/2. ,origin_xd+w_l/2.,num).tolist()
    x_w =  np.linspace(origin_xd, origin_xd,num+1).tolist()

    y_fu = np.linspace(origin_yd + h - t_f_u/2., origin_yd + h - t_f_u/2.,num).tolist()
    y_fl = np.linspace(origin_yd + t_f_l/2., origin_yd + t_f_l/2.,num).tolist()
    y_w =  np.linspace(origin_yd + t_f_l/2.,origin_yd + h - t_f_u/2.,num+1).tolist()

    #Define shear flow discretization by equivalent area determination and flow calculation
    # sigma_node = [(y_w[x] - y_w[0] - (y_w[-1] - y_cdg))/(y_w[-1] - y_cdg) for x in range(len(y_w))]
    # print sigma_node


    # Br_fu = [2*t_f_u*w_u/(num - 1)/2. + t_w*(y_w[-1] - y_w[-2])/6.*(2 + sigma_node[-2]/sigma_node[-1]) if x is x_fu[len(x_fu)/2] else t_f_u*w_u/(num - 1)/2. if x is x_fu[0] else t_f_u*w_u/(num - 1)/2. if x is x_fu[-1] else 2*t_f_u*w_u/(num - 1)/2. for x in x_fu]
    # print Br_fu
    # Br_fl = [2*t_f_l*w_l/(num - 1)/2. + t_w*(y_w[1] -  y_w[0])/6.*(2 +  sigma_node[1]/sigma_node[0])   if x is x_fl[len(x_fl)/2] else t_f_l*w_l/(num - 1)/2. if x is x_fl[0] else t_f_l*w_l/(num - 1)/2. if x is x_fl[-1] else 2*t_f_l*w_l/(num - 1)/2. for x in x_fl]
    # print Br_fl
    # Br_w = [t_w*(y_w[x+1] - y_w[x])/6.*(2 + sigma_node[x]/sigma_node[x+1]) + t_w*(y_w[x+2] - y_w[x+1])/6.*(2 + sigma_node[x+2]/sigma_node[x+1]) for x in range(len(y_w) - 2)]
    # Br_w.insert(0,Br_fl[len(x_fl)/2])
    # Br_w.extend([Br_fu[len(x_fu)/2]])
    # print Br_w 


    # q_dis_fu = [- (V_y*I_y - V_x*I_xy)/(I_x*I_y - I_xy**2)*(Br_fu[x]*(y_fu[x] - y_cdg) + Br_fu[x+1]*(y_fu[x+1] - y_cdg)) - (V_x*I_x - V_y*I_xy)/(I_x*I_y - I_xy**2)*(Br_fu[x]*(x_fu[x] - origin_xd) + Br_fu[x+1]*(x_fu[x+1] - origin_xd)) for x in range(len(y_fu) - 1)]
    # q_dis_fl = [- (V_y*I_y - V_x*I_xy)/(I_x*I_y - I_xy**2)*(Br_fl[x]*(y_fl[x] - y_cdg) + Br_fl[x+1]*(y_fl[x+1] - y_cdg)) - (V_x*I_x - V_y*I_xy)/(I_x*I_y - I_xy**2)*(Br_fl[x]*(x_fl[x] - origin_xd) + Br_fl[x+1]*(x_fl[x+1] - origin_xd)) for x in range(len(y_fl) - 1)]
    # q_dis_w = [- (V_y*I_y - V_x*I_xy)/(I_x*I_y - I_xy**2)*(Br_w[x]*(y_w[x] - y_cdg) + Br_w[x+1]*(y_w[x+1] - y_cdg)) - (V_x*I_x - V_y*I_xy)/(I_x*I_y - I_xy**2)*(Br_w[x]*(x_w[x] - origin_xd) + Br_w[x+1]*(x_w[x+1] - origin_xd)) for x in range(len(y_w) - 1)]

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    x_fl = list(itertools.chain.from_iterable(itertools.repeat(x_fl[x], 2) for x in range(len(x_fl))))
    x_w =list(itertools.chain.from_iterable(itertools.repeat(x_w[x], 2) for x in range(len(x_w))))

    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))
    y_fl = list(itertools.chain.from_iterable(itertools.repeat(y_fl[x], 2) for x in range(len(y_fl))))
    y_w =list(itertools.chain.from_iterable(itertools.repeat(y_w[x], 2) for x in range(len(y_w))))

    #Append constants at beginning and end of the list to offset shape as much as the shear flow value at the intersection with the following element
    y_flange_u = np.linspace(0.,q_flange_u,num/2 + 1).tolist() + np.linspace(q_flange_u,0.,num/2 + 1).tolist()[1:]
    y_flange_l = np.linspace(0.,q_flange_l,num/2 + 1).tolist() + np.linspace(q_flange_l,0.,num/2 + 1).tolist()[1:]

    x_web = lambda y: a_flow*(y)**2 + b_flow*y + c_flow

    #Define discretized shear flow by averaging calculated flow on non-idealized cross section in flanges and web
    q_fu_s = [(y_flange_u[x+1] + y_flange_u[x])/2 for x in range(len(y_flange_u) - 1)]
    q_fl_s = [(y_flange_l[x+1] + y_flange_l[x])/2 for x in range(len(y_flange_l) - 1)]
    q_w_s =  [integrate.quad(x_web,y_w[2*x],y_w[2*x+2])[0]/(y_w[2*x + 2] - y_w[2*x]) for x in range(len(y_w)/2 -1)]

    q_fu = [0.] + [q_fu_s[x/2] for x in range(len(x_fu)-2)] + [0.]
    q_fl = [0.] + [q_fl_s[x/2] for x in range(len(x_fl)-2)] + [0.]
    q_w = [0.] + [q_w_s[x/2] for x in range(len(x_w)-2)] + [0.]

    q_fu_pl = [x + origin_y + h + h/14.*x/max([q_flange_u,q_flange_l]) for x in q_fu]
    q_fl_pl = [x + origin_y  - h/14.*x/max([q_flange_u,q_flange_l]) for x in q_fl]
    q_w_pl = [origin_xd -max([w_u/2.,w_l/2.])- h/14.*x/max([q_flange_u,q_flange_l]) for x in q_w]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]
    q_fl[0] = q_fl_s[0]
    q_fl[-1] = q_fl_s[-1]
    q_w[0] = q_w_s[0]
    q_w[-1] = q_w_s[-1]

    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x = x_fl + x_w + x_fu,
        y = y_fl + y_w + y_fu,
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'none',
        line=dict(color='grey'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the upper flange
    trace2b = go.Scatter(
        x= x_fu,
        y= q_fu_pl,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_fu],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the lower flange
    trace3b = go.Scatter(
        x= x_fl,
        y= q_fl_pl,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_fl],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the web
    trace4b = go.Scatter(
        x=q_w_pl,
        y=y_w,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_w],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Add all the plotting entities to a data list that maps the plotting
    data_temp = [trace1, trace1_sc, trace2, trace3, trace4, trace1b, trace2b, trace3b, trace4b]

    data = data + data_temp
    c = c +1

    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = [dict(x = origin_x - w_u/2. + w_u/10.,y = origin_y + h - t_f_u/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-w_u/4,  ay=0),
                       dict(x = origin_x - w_u/2. + w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-2.5*w_u/6., ay=0),
                       dict(x = origin_x - w_u/2. + w_u/2. - t_w/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-0.7*w_u, ay=0),
                       dict(x = origin_x - w_u/2. +w_u - w_u/10., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=w_u/4, ay=0),
                       dict(x = origin_x - w_u/2. + w_u - w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=2.5*w_u/6., ay=0),
                       dict(x = origin_x - w_u/2. + w_u/2. + t_w/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0.7*w_u, ay=0),
                       
                       dict(x = origin_x - w_l/2. + w_l/10.,y = origin_y + t_f_l/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-w_l/4,  ay=0),
                       dict(x = origin_x - w_l/2. + w_l/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-2.5*w_l/6., ay=0),
                       dict(x = origin_x - w_l/2. + w_l/2. - t_w/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-0.7*w_l, ay=0),
                       dict(x = origin_x - w_l/2. +w_l - w_l/10., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=w_l/4, ay=0),
                       dict(x = origin_x - w_l/2. + w_l - w_l/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=2.5*w_l/6., ay=0),
                       dict(x = origin_x - w_l/2. + w_l/2. + t_w/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0.7*w_l, ay=0),

                       dict(x = origin_x, y = origin_y + t_f_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -0.55*h ),
                       dict(x = origin_x, y = origin_y + 0.775*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -0.55*h ),
                       dict(x = origin_x, y = origin_y + 0.225*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -0.65*h ),
                       dict(x = origin_x, y = origin_y + 0.6*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -0.65*h ),
                       dict(x = origin_x, y = origin_y + 0.4*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -0.75*h ),

                       dict(x = origin_x, y = y_cdg + t_w/4., xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 16 ),showarrow = True, arrowhead=2, arrowsize=2, arrowwidth=1, arrowcolor='black', ax=0., ay= -h/2., xanchor = "center" ) ],
        #Define axes visibility properties
        xaxis=dict(
            autorange=False,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            range=[-1.5*max([w_u,w_l]), 1.5*max([w_u,w_l])],
            showticklabels=False
        ),
        yaxis=dict(
            autorange=False,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            range=[-0.2*h, 1.2*h],
            scaleanchor="x", 
            scaleratio=1,
            showticklabels=False
        )
    )

    return (data, layout)

opts = [{'label': 'I Shape', 'value': 'I'},]

data, layout = gen_data_i_beam(V_y,V_x,t_w,t_f_u,t_f_l,h,w_u,w_l,num)

fig = go.Figure(data= data, layout = layout)

app.layout = html.Div([
    html.Div([html.Div([
        html.Label("Choose cross-section type"),
        dcc.Dropdown(id = 'opt', options = opts,
                    value = 'I')
            ], style = {'width': '168px',
                        'fontSize' : '12px',
                        'padding-left' : '100px',
                        'display': 'table-cell'}),
    html.Div([
    html.Label("Define shear load (V) [N]"),
    dcc.Input(
        id='inp_V',
        placeholder='Enter a value...',
        type='number',
        value=V_y)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'}),
    ]),
    html.Div([
    html.Label("Define section height (h) [mm]"),
    dcc.Input(
        id='inp_h',
        placeholder='Enter a value...',
        type='number',
        value=h)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '100px',
                    'display': 'table-cell'}),
    html.Div([
    html.Label("Define top width (w1) [mm]"),
    dcc.Input(
        id='inp_w1',
        placeholder='Enter a value...',
        type='number',
        value=w_u)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'}),
    html.Div([
    html.Label("Define bottom width (w2) [mm]"),
    dcc.Input(
        id='inp_w2',
        placeholder='Enter a value...',
        type='number',
        value=w_l)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'}),
    html.Div([
    html.Label("Define top thickness (t1) [mm]"),
    dcc.Input(
        id='inp_t1',
        placeholder='Enter a value...',
        type='number',
        value=t_f_u)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'}),
    html.Div([
    html.Label("Define bottom thick. (t2) [mm]"),
    dcc.Input(
        id='inp_t2',
        placeholder='Enter a value...',
        type='number',
        value=t_f_l)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'}),
    html.Div([
    html.Label("Define web thickness (t3) [mm]"),
    dcc.Input(
        id='inp_t3',
        placeholder='Enter a value...',
        type='number',
        value=t_w)
        ], style = {'width': '50px',
                    'fontSize' : '12px',
                    'padding-left' : '20px',
                    'display': 'table-cell'}),
    dcc.Graph(
        id='flow',
        figure=fig
    ),
    html.Div([
    html.Label('Define discretization'),
    dcc.Slider(id='slider',min=3,max=21,step=2,value=7,
    marks={i: '{}'.format(i) if i == 3 else str(i) for i in np.arange(3, 23,2)},
    ),],
        style = {'width' : '80%',
                'fontSize' : '12px',
                'padding-left' : '100px',
                'display': 'inline-block'}
        
        ),
        html.Div(id='slider-output-container')
    ])


@app.callback(
    dash.dependencies.Output('flow', 'figure'),
    [dash.dependencies.Input('slider', 'value'),
    dash.dependencies.Input('opt', 'value'),
    dash.dependencies.Input('inp_V', 'value'),
    dash.dependencies.Input('inp_h', 'value'),
    dash.dependencies.Input('inp_w1', 'value'),
    dash.dependencies.Input('inp_w2', 'value'),
    dash.dependencies.Input('inp_t1', 'value'),
    dash.dependencies.Input('inp_t2', 'value'),
    dash.dependencies.Input('inp_t3', 'value')])

def update_output(input1,input2,input3,input4,input5,input6,input7,input8,input9):
    data, layout = gen_data_i_beam(input3,V_x,input9,input7,input8,input4,input5,input6,input1)
    fig = go.Figure(data = data, layout = layout)
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
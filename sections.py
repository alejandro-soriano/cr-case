import numpy as np
import math
import scipy.optimize as opt
import scipy.integrate as integrate
import scipy.linalg as linalg
import itertools

import plotly.graph_objs as go

#Functions that determines and organizes the data used for plotting

#Define equation of a parabola for shear flow field plotting
def eq_parab(p,e1,e2,vert): #Function used to calculate the parabola constants
    a, b, c = p
    eq1 = -e1[0] + a*e1[1]**2 + b*e1[1] + c
    # eq2 = -e2[0] + a*e2[1]**2 + b*e2[1] + c
    eq2 = -vert[0] + a*vert[1]**2 + b*vert[1] + c
    eq3 = 2*a*vert[1] + b
    return (eq1,eq2,eq3)


def gen_data_i_beam(V_y,V_x,t_w,t_f_u,t_f_l,h,w_u,w_l,num):

    #Calculate geometrical properties of section
    A_tot = w_u*t_f_u + t_w*(h - t_f_u - t_f_l) + w_l*t_f_l #mm2
    y_cdg = (w_u*t_f_u*(h - t_f_u/2.) + t_w*(h - t_f_u - t_f_l)*h/2. + w_l*t_f_l**2/2)/A_tot #mm w.r.t. bottom line
    I_x = 1./12.*w_u*t_f_u**3 + 1./12.*w_l*t_f_l**3 + w_u*t_f_u*(h - y_cdg - t_f_u/2.)**2 + w_l*t_f_l*(y_cdg - t_f_l/2.)**2 + 1./12.*t_w*(h - t_f_u - t_f_l)**3 +(h - t_f_u - t_f_l)*t_w*(h/2. - y_cdg)**2  #mm4
    I_y = 1./12.*t_f_u*w_u**3 + 1./12.*t_f_l*w_l**3 + 1./12.*(h-t_f_u-t_f_l)*t_w**3
    I_xy = 0.


    #Compute shear flows on section
    q_flange_u = V_y*(t_f_u*w_u/2.*(h - y_cdg - t_f_u/2.))/I_x #N/mm
    q_flange_l = V_y*(t_f_l*w_l/2.*(y_cdg - t_f_l/2.))/I_x #N/mm


    Q_w_u = t_f_u*w_u*(h - y_cdg - t_f_u/2.) + (h - y_cdg - t_f_u)**2/2*t_w #mm3
    Q_w_l = t_f_l*w_l*(y_cdg - t_f_l/2.) + (y_cdg - t_f_l)**2/2*t_w #mm3

    q_web = V_y*Q_w_u/I_x #N/mm @ cog

    q_web_fun = lambda y: V_y*(t_f_u*w_u*(h - y_cdg - t_f_u/2.) + (h - y - t_f_u)*t_w*(((h - t_f_u) + y)/2. - y_cdg))/I_x

    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    if max([w_u,w_l]) > h:
        origin_x = -(h/max([w_u,w_l])*(h - 0.75*max([w_u,w_l])) + 0.75*max([w_u,w_l]))
    else:
        origin_x = -1.*h
    origin_y = 0.

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1 = go.Scatter(
        x=[origin_x - w_l/2., origin_x + w_l/2. ,origin_x+w_l/2., origin_x + t_w/2.,origin_x + t_w/2.,origin_x+w_u/2.,origin_x+w_u/2.,origin_x - w_u/2.,origin_x - w_u/2., origin_x-t_w/2.,origin_x-t_w/2.,origin_x - w_l/2.,origin_x - w_l/2.],
        y=[origin_y,origin_y,origin_y+t_f_l,origin_y+t_f_l,origin_y+h-t_f_u,origin_y+h-t_f_u,origin_y+h,origin_y+h,origin_y+h-t_f_u,origin_y+h-t_f_u,origin_y+t_f_l,origin_y+t_f_l,origin_y],
        mode='lines',
        name=r'I Beam <br> V = {0:.2f} N <br> h = {1:.0f} mm <br> w1 = {2:.0f} mm <br> w2 = {3:.0f} mm <br> t1 = {4:.0f} mm <br> t2 = {5:.0f} mm <br> t3 = {6:.0f} mm <br> A = {7:.0f} mm2 <br> I_x = {8:.2f} mm4 <br> I_y = {9:.0f} mm4'.format(V_y,h,w_u,w_l,t_f_u,t_f_l,t_w,A_tot,I_x,I_y),
        hoverinfo='text',
        hoveron='fills',
        fill='toself', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    x_sc = 0.
    y_sc = h - h*t_f_l/I_y*w_l**3/12

    #Entity that plots the location of the shear load application (shear center)
    trace1_sc = go.Scatter(
        x = [origin_x],
        y = [origin_y + y_sc - t_f_l/2.*(y_sc - h/2.)/h*2.],
        name = 'Shear Center',
        hovertext = r'Shear Center <br> x_sc = {0:.2f} <br> y_sc = {1:.2f} <br> w.r.t.bottom & plane of symmetry'.format(x_sc,round(y_sc - t_f_l/2.*(y_sc - h/2.)/h*2.,1)),
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'blue', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the location of the center of gravity
    trace1_cg = go.Scatter(
        x = [origin_x],
        y = [origin_y + y_cdg],
        hovertext = r'Center of Gravity <br> x_cog = {0:.2f} <br> y_cog = {1:.2f} <br> w.r.t. bottom & plane of symmetry'.format(0.,y_cdg),
        name = 'Center of Gravity',
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'purple', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the neutral axis
    trace1_na = go.Scatter(
        x = np.linspace(origin_x - max([w_u,w_l])/4.,origin_x+ max([w_u,w_l])/4.,21),
        y = np.linspace(origin_y + y_cdg, origin_y + y_cdg, 21),
        name = 'Neutral Axis',
        mode='lines',
        hoveron = 'points',
        hoverinfo = 'name',
        line = dict(
            color = ('purple'),
            dash = 'dashdot'))

    #Determination of the shear flow fields in the upper and lower flanges of the I-section
    y_flange_u = np.linspace(0.,q_flange_u,11).tolist() + np.linspace(q_flange_u,0.,11).tolist()
    y_flange_l = np.linspace(0.,q_flange_l,11).tolist() + np.linspace(q_flange_l,0.,11).tolist()

    #Entity that plots the upper flange shear flow
    trace2 = go.Scatter(
        x= np.linspace(origin_x - w_u/2.,origin_x,11).tolist() + np.linspace(origin_x,origin_x + w_u/2.,11).tolist(),
        y= np.linspace(origin_y + h,origin_y + h + h/5.*q_flange_u/max([q_flange_u,q_flange_l]),11).tolist() + np.linspace(origin_y + h + h/5.*q_flange_u/max([q_flange_u,q_flange_l]),origin_y + h,11).tolist(),
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
        y= np.linspace(origin_y ,origin_y - h/5.*q_flange_l/max([q_flange_u,q_flange_l]),11).tolist() + np.linspace(origin_y - h/5.*q_flange_l/max([q_flange_u,q_flange_l]),origin_y ,11).tolist(),
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
    #Points through which the parabola-shaped shear flow in the web goes in order to define the appropriate tags
    e1_flow = [q_flange_l,t_f_l/2.]
    e2_flow = [q_flange_u,h - t_f_u/2.]
    vert_flow = [q_web, y_cdg]

    #Solve for the parabola constants
    a_flow, b_flow, c_flow = opt.fsolve(eq_parab,[1/180.,-5./18.,713./36.],args=(e1_flow,e2_flow,vert_flow))

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow = np.linspace(origin_y+t_f_l/2.,origin_y+h-t_f_u/2.,21)
    x_par_flow = q_web_fun(y_par_flow)

    x_par_coord = [origin_x - max([w_u/2.,w_l/2.]) - x_par_flow[x]*h/5./max([q_flange_u,q_flange_l]) for x in range(len(x_par_flow))]
    x_par_flow = x_par_flow.tolist()
    x_par_flow.append(q_flange_u)
    x_par_flow.insert(0,q_flange_l)

    #Points through which the parabola-shaped shear flow in the web goes
    e1_coord = [origin_x - max([w_u/2.,w_l/2.]) - h/5.*q_flange_l/max([q_flange_u,q_flange_l]),t_f_l/2.]
    e2_coord = [origin_x - max([w_u/2.,w_l/2.])- h/5.*q_flange_u/max([q_flange_u,q_flange_l]),h - t_f_u/2.]
    vert_coord = [origin_x -max([w_u/2.,w_l/2.])- h/5.*q_web/max([q_flange_u,q_flange_l]),y_cdg]

    #Solve for the parabola constants
    a_coord, b_coord, c_coord = opt.fsolve(eq_parab,[1/180.,-5./18.,713./36.],args=(e1_coord,e2_coord,vert_coord))

    y_par = np.linspace(origin_y+t_f_l/2.,origin_y+h-t_f_u/2.,21)

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par = y_par.tolist()
    y_par.append(origin_y+h-t_f_u/2.)
    y_par.insert(0,origin_y+t_f_l/2.)

    x_par_coord.append(origin_x - max([w_u/2.,w_l/2.]))
    x_par_coord.insert(0,origin_x - max([w_u/2.,w_l/2.]))

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
    if max([w_u,w_l]) > h:
        origin_xd = (h/max([w_u,w_l])*(h - 0.75*max([w_u,w_l])) + 0.75*max([w_u,w_l]))
    else:
        origin_xd = 1.*h
    origin_yd = 0.

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = np.linspace(origin_xd - w_u/2. ,origin_xd+w_u/2.,num).tolist()
    x_fl = np.linspace(origin_xd - w_l/2. ,origin_xd+w_l/2.,num).tolist()
    x_w =  np.linspace(origin_xd, origin_xd,num+1).tolist()

    y_fu = np.linspace(origin_yd + h - t_f_u/2., origin_yd + h - t_f_u/2.,num).tolist()
    y_fl = np.linspace(origin_yd + t_f_l/2., origin_yd + t_f_l/2.,num).tolist()
    y_w =  np.linspace(origin_yd + t_f_l/2.,origin_yd + h - t_f_u/2.,num+1).tolist()

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    x_fl = list(itertools.chain.from_iterable(itertools.repeat(x_fl[x], 2) for x in range(len(x_fl))))
    x_w =list(itertools.chain.from_iterable(itertools.repeat(x_w[x], 2) for x in range(len(x_w))))

    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))
    y_fl = list(itertools.chain.from_iterable(itertools.repeat(y_fl[x], 2) for x in range(len(y_fl))))
    y_w =list(itertools.chain.from_iterable(itertools.repeat(y_w[x], 2) for x in range(len(y_w))))

    #Append constants at beginning and end of the list to offset shape as much as the shear flow value at the intersection with the following element
    y_flange_u = np.linspace(0.,q_flange_u,int(num/2 + 1)).tolist() + np.linspace(q_flange_u,0.,int(num/2 + 1)).tolist()[1:]
    y_flange_l = np.linspace(0.,q_flange_l,int(num/2 + 1)).tolist() + np.linspace(q_flange_l,0.,int(num/2 + 1)).tolist()[1:]

    x_web = q_web_fun

    #Define discretized shear flow by averaging calculated flow on non-idealized cross section in flanges and web
    q_fu_s = [(y_flange_u[x+1] + y_flange_u[x])/2 for x in range(len(y_flange_u) - 1)]
    q_fl_s = [(y_flange_l[x+1] + y_flange_l[x])/2 for x in range(len(y_flange_l) - 1)]
    q_w_s =  [integrate.quad(x_web,y_w[2*x],y_w[2*x+2])[0]/(y_w[2*x + 2] - y_w[2*x]) for x in range(int(len(y_w)/2 -1))]

    q_w_total = integrate.quad(x_web,y_w[0],y_w[-1])[0]

    q_fu = [0.] + [q_fu_s[int(x/2)] for x in range(int(len(x_fu)-2))] + [0.]
    q_fl = [0.] + [q_fl_s[int(x/2)] for x in range(int(len(x_fl)-2))] + [0.]
    q_w = [0.] + [q_w_s[int(x/2)] for x in range(int(len(x_w)-2))] + [0.]

    q_fu_pl = [x + origin_y + h + h/5.*x/max([q_flange_u,q_flange_l]) for x in q_fu]
    q_fl_pl = [x + origin_y  - h/5.*x/max([q_flange_u,q_flange_l]) for x in q_fl]
    q_w_pl = [origin_xd -max([w_u/2.,w_l/2.])- h/5.*x/max([q_flange_u,q_flange_l]) for x in q_w]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]
    q_fl[0] = q_fl_s[0]
    q_fl[-1] = q_fl_s[-1]
    q_w[0] = q_w_s[0]
    q_w[-1] = q_w_s[-1]

    #Calculate the corresponding boom areas based on the ratios of shear flow between two adjacent booms
    A_fu = np.zeros(num)
    A_fl = np.zeros(num)
    A_w = np.zeros(num+1)

    y_w_A =  np.linspace(origin_yd,origin_yd + h - t_f_u/2.,num+1).tolist()
    x_fu_A = np.linspace(origin_xd ,origin_xd+w_u/2.,num).tolist()
    x_fl_A = np.linspace(origin_xd ,origin_xd+w_l/2.,num).tolist()

    q_fupper = lambda s: q_flange_u/w_u*2*s

    for x in range(0,math.ceil(num/2) - 1):
        A_fu[x] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x] - origin_xd)]))
        A_fu[x +1] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x +1] - origin_xd)]))
        
        A_fu[num - x -1] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x] - origin_xd)]))
        A_fu[num - x -2] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x +1] - origin_xd)]))
        
        A_fl[x] += t_f_l*w_l/(num -1)/6.*(2 + q_fupper(x_fl_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x] - origin_xd)]))
        A_fl[x +1] += t_f_l*w_l/(num -1)/6.*(2 + q_fupper(x_fl_A[x] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x +1] - origin_xd)]))
        
        A_fl[num - x -1] += t_f_l*w_l/(num -1)/6.*(2 + q_fupper(x_fl_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x] - origin_xd)]))
        A_fl[num - x -2] += t_f_l*w_l/(num -1)/6.*(2 + q_fupper(x_fl_A[x] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x +1] - origin_xd)]))


    for x in range(0, num):
        A_w[x] += t_w*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x+1] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x] - origin_yd)]))
        A_w[x+1] += t_w*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x+1] - origin_yd)]))

      
    A_fu[math.ceil(num/2) - 1] = A_fu[math.ceil(num/2) - 1] + A_w[-1]
    A_fl[math.ceil(num/2) - 1] = A_fl[math.ceil(num/2) - 1] + A_w[0]
    A_w[-1] = A_fu[math.ceil(num/2) - 1]
    A_w[0] = A_fl[math.ceil(num/2) - 1]

    A_fu = A_fu.tolist()
    A_fl = A_fl.tolist()
    A_w = A_w.tolist()

    res_A = A_fl + A_w + A_fu
    res_A = list(itertools.chain.from_iterable(itertools.repeat(res_A[x], 2) for x in range(len(res_A))))

    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x = x_fl + x_w + x_fu,
        y = y_fl + y_w + y_fu,
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['Ab = {0:.0f} mm2'.format(k) for k in res_A],
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
    data_temp = [trace1, trace1_sc, trace1_cg,trace1_na, trace2, trace3, trace4, trace1b, trace2b, trace3b, trace4b]

    data = data + data_temp
    c = c +1

    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = [dict(x = origin_x - w_u/2. + w_u/10.,y = origin_y + h - t_f_u/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-18,  ay=0),
                       dict(x = origin_x - w_u/2. + w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-30, ay=0),
                       dict(x = origin_x - w_u/2. + w_u/2. - t_w/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-45, ay=0),
                       dict(x = origin_x - w_u/2. +w_u - w_u/10., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=18, ay=0),
                       dict(x = origin_x - w_u/2. + w_u - w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=30, ay=0),
                       dict(x = origin_x - w_u/2. + w_u/2. + t_w/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=45, ay=0),
                       
                       dict(x = origin_x - w_l/2. + w_l/10.,y = origin_y + t_f_l/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-18,  ay=0),
                       dict(x = origin_x - w_l/2. + w_l/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-30, ay=0),
                       dict(x = origin_x - w_l/2. + w_l/2. - t_w/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-45, ay=0),
                       dict(x = origin_x - w_l/2. +w_l - w_l/10., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=18, ay=0),
                       dict(x = origin_x - w_l/2. + w_l - w_l/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=30, ay=0),
                       dict(x = origin_x - w_l/2. + w_l/2. + t_w/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=45, ay=0),

                       dict(x = origin_x, y = origin_y + t_f_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x, y = origin_y + 0.85*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x, y = origin_y + 0.2*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x, y = origin_y + 0.65*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x, y = origin_y + 0.37*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),

                       dict(x = origin_x, y = origin_y + y_sc - t_f_l/2.*(y_sc - h/2.)/h*2. + 0.05*h, xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 14 ),showarrow = True, arrowhead=1, arrowsize=1.5, arrowwidth=1, arrowcolor='black', ax=0., ay= -40, xanchor = "center" ) ],
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

def gen_data_t_beam(V_y,V_x,t_w,t_f_u,h,w_u,num):

    #Calculate geometrical properties of section
    A_tot = w_u*t_f_u + t_w*(h - t_f_u)  #mm2
    y_cdg = (w_u*t_f_u*(h - t_f_u/2.) + t_w*(h - t_f_u)**2/2.)/A_tot #mm w.r.t. bottom line
    x_cdg = 0.
    I_x = 1./12.*w_u*t_f_u**3 + w_u*t_f_u*(h - y_cdg - t_f_u/2.)**2 + 1./12.*t_w*(h - t_f_u)**3 +(h - t_f_u)*t_w*(h/2. - y_cdg)**2  #mm4
    I_y = 1./12.*t_f_u*w_u**3 + 1./12.*(h-t_f_u)*t_w**3
    I_xy = 0.

    #Compute maximum shear flows on section
    q_flange_u = V_y*(t_f_u*w_u/2.*(h - y_cdg - t_f_u/2.))/I_x #N/mm

    Q_w_u = t_f_u*w_u*(h - y_cdg - t_f_u/2.) + (h - y_cdg - t_f_u)**2/2*t_w #mm3

    q_web = V_y*Q_w_u/I_x #N/mm @ cog

    q_web_fun = lambda y: V_y*(t_f_u*w_u*(h - y_cdg - t_f_u/2.) + (h - y - t_f_u)*t_w*(((h - t_f_u) + y)/2. - y_cdg))/I_x

    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    if max([w_u]) > h:
        origin_x = -(h/w_u*(h - 0.75*w_u) + 0.75*w_u)
    else:
        origin_x = -1.*h
    origin_y = 0.

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1 = go.Scatter(
        x=[origin_x + t_w/2.,origin_x + t_w/2.,origin_x+w_u/2.,origin_x+w_u/2.,origin_x - w_u/2.,origin_x - w_u/2., origin_x-t_w/2.,origin_x-t_w/2.,origin_x + t_w/2.],
        y=[origin_y,origin_y+h-t_f_u,origin_y+h-t_f_u,origin_y+h,origin_y+h,origin_y+h-t_f_u,origin_y+h-t_f_u,origin_y,origin_y],
        mode='lines',
        name=r'T Beam <br> V = {0:.2f} N <br> h = {1:.0f} mm <br> w1 = {2:.0f} mm <br> t1 = {4:.0f} mm <br> t3 = {6:.0f} mm <br> A = {7:.0f} mm2 <br> I_x = {8:.2f} mm4 <br> I_y = {9:.0f} mm4'.format(V_y,h,w_u,w_u,t_f_u,t_f_u,t_w,A_tot,I_x,I_y),
        hoverinfo='text',
        hoveron='fills',
        fill='toself', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    x_sc = 0.
    y_sc = h - t_f_u/2.

    #Entity that plots the location of the shear load application (shear center)
    trace1_sc = go.Scatter(
        x = [origin_x],
        y = [origin_y + h - t_f_u/2.],
        name = 'Shear Center',
        hovertext = r'Shear Center <br> x_sc = {0:.2f} <br> y_sc = {1:.2f} <br> w.r.t.bottom & plane of symmetry'.format(x_sc,y_sc),
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'blue', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the location of the center of gravity
    trace1_cg = go.Scatter(
        x = [origin_x],
        y = [origin_y + y_cdg],
        hovertext = r'Center of Gravity <br> x_cog = {0:.2f} <br> y_cog = {1:.2f} <br> w.r.t. bottom & plane of symmetry'.format(0.,y_cdg),
        name = 'Center of Gravity',
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'purple', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the neutral axis
    trace1_na = go.Scatter(
        x = np.linspace(origin_x - w_u/4.,origin_x+ w_u/4.,21),
        y = np.linspace(origin_y + y_cdg, origin_y + y_cdg, 21),
        name = 'Neutral Axis',
        mode='lines',
        hoveron = 'points',
        hoverinfo = 'name',
        line = dict(
            color = ('purple'),
            dash = 'dashdot'))


    #Determination of the shear flow fields in the upper and lower flanges of the I-section
    y_flange_u = np.linspace(0.,q_flange_u,11).tolist() + np.linspace(q_flange_u,0.,11).tolist()

    #Entity that plots the upper flange shear flow
    trace2 = go.Scatter(
        x= np.linspace(origin_x - w_u/2.,origin_x,11).tolist() + np.linspace(origin_x,origin_x + w_u/2.,11).tolist(),
        y= np.linspace(origin_y + h,origin_y + h + h/5.*q_flange_u/max([q_flange_u]),11).tolist() + np.linspace(origin_y + h + h/5.*q_flange_u/max([q_flange_u]),origin_y + h,11).tolist(),
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in y_flange_u],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Computation of parabola constats for plotting the shear flow on web
    #Points through which the parabola-shaped shear flow in the web goes in order to define the appropriate tags
    e1_flow = [0.,origin_y]
    e2_flow = [q_flange_u,origin_y + h - t_f_u/2.]
    vert_flow = [q_web, y_cdg]

    #Solve for the parabola constants
    a_flow, b_flow, c_flow = opt.fsolve(eq_parab,[1/180.,-5./18.,713./36.],args=(e1_flow,e2_flow,vert_flow))

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow = np.linspace(origin_y,origin_y+h-t_f_u/2.,21)
    x_par_flow = q_web_fun(y_par_flow)

    x_par_coord = [origin_x - max([w_u/2.,w_u/2.]) - x_par_flow[x]*h/5./max([q_flange_u]) for x in range(len(x_par_flow))]

    x_par_flow = x_par_flow.tolist()
    x_par_flow.append(q_flange_u)
    x_par_flow.insert(0,0.)

    #Points through which the parabola-shaped shear flow in the web goes
    e1_coord = [origin_x - max([w_u/2.,w_u/2.]) - h/5.*0./max([q_flange_u]),origin_y]
    e2_coord = [origin_x - max([w_u/2.,w_u/2.])- h/5.*q_flange_u/max([q_flange_u]),origin_y + h - t_f_u/2.]
    vert_coord = [origin_x -max([w_u/2.,w_u/2.])- h/5.*q_web/max([q_flange_u]),origin_y + y_cdg]
    
    #Solve for the parabola constants
    a_coord, b_coord, c_coord = opt.fsolve(eq_parab,[1/180.,-5./18.,713./36.],args=(e1_coord,e2_coord,vert_coord))
    y_par = np.linspace(origin_y,origin_y+h-t_f_u/2.,21)
    x_par_coord = a_coord*(y_par)**2 + b_coord*y_par + c_coord

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par = y_par.tolist()
    y_par.append(origin_y+h-t_f_u/2.)
    y_par.insert(0,origin_y)

    x_par_coord = x_par_coord.tolist()
    x_par_coord.append(origin_x - w_u/2.)
    x_par_coord.insert(0,origin_x -  w_u/2.)

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
    if max([w_u]) > h:
        origin_xd = (h/max([w_u,w_u])*(h - 0.75*max([w_u,w_u])) + 0.75*max([w_u,w_u]))
    else:
        origin_xd = 1.*h

    origin_yd = 0.

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = np.linspace(origin_xd - w_u/2. ,origin_xd+w_u/2.,num).tolist()
    x_w =  np.linspace(origin_xd, origin_xd,num+1).tolist()

    y_fu = np.linspace(origin_yd + h - t_f_u/2., origin_yd + h - t_f_u/2.,num).tolist()
    y_w =  np.linspace(origin_yd,origin_yd + h - t_f_u/2.,num+1).tolist()

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    x_w =list(itertools.chain.from_iterable(itertools.repeat(x_w[x], 2) for x in range(len(x_w))))

    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))
    y_w =list(itertools.chain.from_iterable(itertools.repeat(y_w[x], 2) for x in range(len(y_w))))

    #Append constants at beginning and end of the list to offset shape as much as the shear flow value at the intersection with the following element
    y_flange_u = np.linspace(0.,q_flange_u,int(num/2 + 1)).tolist() + np.linspace(q_flange_u,0.,int(num/2 + 1)).tolist()[1:]

    x_web = lambda y: a_flow*(y)**2 + b_flow*y + c_flow

    #Define discretized shear flow by averaging calculated flow on non-idealized cross section in flanges and web
    q_fu_s = [(y_flange_u[x+1] + y_flange_u[x])/2 for x in range(len(y_flange_u) - 1)]
    q_w_s =  [integrate.quad(x_web,y_w[2*x],y_w[2*x+2])[0]/(y_w[2*x + 2] - y_w[2*x]) for x in range(int(len(y_w)/2 -1))]

    q_fu = [0.] + [q_fu_s[int(x/2)] for x in range(int(len(x_fu)-2))] + [0.]
    q_w = [0.] + [q_w_s[int(x/2)] for x in range(int(len(x_w)-2))] + [0.]

    q_fu_pl = [x + origin_y + h + h/5.*x/max([q_flange_u]) for x in q_fu]
    q_w_pl = [origin_xd -max([w_u/2.,w_u/2.])- h/5.*x/max([q_flange_u]) for x in q_w]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]
    q_w[0] = q_w_s[0]
    q_w[-1] = q_w_s[-1]

    #Calculate the corresponding boom areas based on the ratios of shear flow between two adjacent booms
    A_fu = np.zeros(num)
    A_w = np.zeros(num+1)

    y_w_A =  np.linspace(origin_yd,origin_yd + h - t_f_u/2.,num+1).tolist()
    x_fu_A = np.linspace(origin_xd ,origin_xd+w_u/2.,num).tolist()

    q_fupper = lambda s: q_flange_u/w_u*2*s

    for x in range(0,math.ceil(num/2) - 1):
        A_fu[x] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x] - origin_xd)]))
        A_fu[x +1] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x +1] - origin_xd)]))
        
        A_fu[num - x -1] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x] - origin_xd)]))
        A_fu[num - x -2] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x +1] - origin_xd)]))

    for x in range(0, num):
        A_w[x] += t_w*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x+1] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x] - origin_yd)]))
        A_w[x+1] += t_w*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x+1] - origin_yd)]))
    
    A_fu[math.ceil(num/2) - 1] = A_fu[math.ceil(num/2) - 1] + A_w[-1]
    A_w[-1] = A_fu[math.ceil(num/2) - 1]

    A_fu = A_fu.tolist()
    A_w = A_w.tolist()

    res_A = A_w + A_fu
    res_A = list(itertools.chain.from_iterable(itertools.repeat(res_A[x], 2) for x in range(len(res_A))))

    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x =  x_w + x_fu,
        y =  y_w + y_fu,
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['Ab = {0:.0f} mm2'.format(k) for k in res_A],
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
    data_temp = [trace1, trace1_sc, trace1_cg, trace1_na, trace2, trace4, trace1b, trace2b, trace4b]

    data = data + data_temp
    c = c +1

    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = [dict(x = origin_x - w_u/2. + w_u/10.,y = origin_y + h - t_f_u/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-18,  ay=0),
                       dict(x = origin_x - w_u/2. + w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-30, ay=0),
                       dict(x = origin_x - w_u/2. + w_u/2. - t_w/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-45, ay=0),
                       dict(x = origin_x - w_u/2. +w_u - w_u/10., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=18, ay=0),
                       dict(x = origin_x - w_u/2. + w_u - w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=30, ay=0),
                       dict(x = origin_x - w_u/2. + w_u/2. + t_w/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=45, ay=0),
                       
                       dict(x = origin_x, y = origin_y, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x, y = origin_y + 0.85*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x, y = origin_y + 0.2*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x, y = origin_y + 0.65*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x, y = origin_y + 0.37*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),

                       dict(x = origin_x, y = origin_y + y_sc + 0.05*h, xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 14 ),showarrow = True, arrowhead=1, arrowsize=1.5, arrowwidth=1, arrowcolor='black', ax=0., ay= -40, xanchor = "center" ) ],
      #Define axes visibility properties
        xaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            # range=[-1.5*max([w_u,w_l]), 1.5*max([w_u,w_l])],
            showticklabels=False
        ),
        yaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            # range=[-0.2*h, 1.2*h],
            scaleanchor="x", 
            scaleratio=1,
            showticklabels=False
        )
    )

    return (data, layout)

def gen_data_c_beam(V_y,V_x,t_f_u,t_w_r,t_w_l,w,h_r,h_l,num):

    #Calculate geometrical properties of section
    A_tot = w*t_f_u + t_w_l*(h_l - t_f_u) + t_w_r*(h_r - t_f_u)
    y_cdg = (t_f_u*w*w/2. + t_w_r*h_r*w)/(t_f_u*w + t_w_r*h_r + t_w_l*h_l)
    x_cdg = (t_w_l*h_l*h_l/2. + t_w_r*h_r*h_r/2.)/(t_f_u*w + t_w_r*h_r + t_w_l*h_l)

    I_x = t_w_l*h_l*y_cdg**2 + (t_f_u*w**3/12) + t_f_u*w*(w/2. - y_cdg)**2 + t_w_r*h_r*(w - y_cdg)**2
    I_y = (t_w_l*h_l**3/12.) + t_w_l*h_l*(h_l/2. - x_cdg)**2 + t_f_u*w*x_cdg**2 + (t_w_r*h_r**3/12.) + t_w_r*h_r*(h_r/2. - x_cdg)**2
    I_xy = t_w_l*h_l*(h_l/2. - x_cdg)*y_cdg + t_f_u*w*(-w/2. + y_cdg)*(-x_cdg) + t_w_r*h_r*(h_r/2. - x_cdg)*(-w + y_cdg)
   
    #Compute maximum shear flows on section
    q43_y = lambda s: V_y*(I_xy/(I_x*I_y - I_xy**2)*t_w_r*((h_r - x_cdg)*s - s**2/2.) + I_y/(I_x*I_y - I_xy**2)*t_w_r*(w - y_cdg)*s)
    q43_x = lambda s: V_x*(I_x/(I_x*I_y - I_xy**2)*t_w_r*((h_r - x_cdg)*s - s**2/2.) + I_xy/(I_x*I_y - I_xy**2)*t_w_r*(w - y_cdg)*s)

    y_sc = integrate.quad(q43_y,0.,h_r)[0]*w/V_y
    x_sc = integrate.quad(q43_x,0.,h_r)[0]*w/V_x


    q12_y = lambda s: V_y*(-I_xy/(I_x*I_y - I_xy**2)*t_w_l*((h_l - x_cdg)*s - s**2/2.) + I_y/(I_x*I_y - I_xy**2)*t_w_l*(y_cdg)*s)
    q12_x = lambda s: V_x*(-I_x/(I_x*I_y - I_xy**2)*t_w_l*((h_l - x_cdg)*s - s**2/2.) + I_xy/(I_x*I_y - I_xy**2)*t_w_l*(y_cdg)*s)

    q23_y = lambda s: V_y*(I_xy/(I_x*I_y - I_xy**2)*t_f_u*(-(x_cdg)*s) - I_y/(I_x*I_y - I_xy**2)*t_f_u*(-(w - y_cdg)*s + s**2/2.)) + q43_y(h_r)
    q23_x = lambda s: V_x*(I_x/(I_x*I_y - I_xy**2)*t_f_u*(-(x_cdg)*s) - I_xy/(I_x*I_y - I_xy**2)*t_f_u*(-(w - y_cdg)*s + s**2/2.)) + q43_x(h_r)

    q_web_fun_l = lambda y: V_l*(t_f_u*(w/2. + x_sc)*(h_l - y_cdg_l - t_f_u/2.) + (h_l - y - t_f_u)*t_w_l*(((h_l - t_f_u) + y)/2. - y_cdg_l))/I_x_l

    q_web_fun_r = lambda y: V_r*(t_f_u*(w/2. - x_sc)*(h_r - y_cdg_r - t_f_u/2.) + (h_r - y - t_f_u)*t_w_r*(((h_r - t_f_u) + y)/2. - y_cdg_r))/I_x_r

    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    if max([h_r,h_l]) > w:
        origin_x = -(w/max([h_r,h_l])*(w - 0.75*max([h_r,h_l])) + 0.75*max([h_r,h_l]))
    else:
        origin_x = -1.*w
    origin_y = 0.

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1 = go.Scatter(
        x=[origin_x, origin_x + h_l ,origin_x+h_l, origin_x  + t_f_u, origin_x + t_f_u,origin_x+h_r,origin_x+h_r,origin_x,origin_x],
        y=[origin_y,origin_y,origin_y+t_w_l,origin_y+t_w_l,origin_y+w-t_w_r,origin_y+w-t_w_r,origin_y+w,origin_y+w,origin_y],
        mode='lines',
        name=r'C Beam <br> V = {0:.2f} N <br> h1 = {1:.0f} mm <br> h2 = {2:.0f} mm <br> w = {3:.0f} mm <br> t2 = {4:.0f} mm <br> t3 = {5:.0f} mm <br> t4 = {6:.0f} mm <br> A = {7:.0f} mm2 <br> I_x = {8:.0f} mm4 <br> I_y = {9:.0f} mm4 <br> I_xy = {10:.0f} mm4'.format(V_y,h_l,h_r,w,t_f_u,t_w_l,t_w_r,A_tot,I_x,I_y,I_xy),
        hoverinfo='text',
        hoveron='fills',
        fill='toself', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    #Entity that plots the location of the shear load application (shear center)
    trace1_sc = go.Scatter(
        x = [origin_x  - y_sc + t_f_u/2.],
        y = [origin_y + x_sc - t_w_r/2.*(x_sc - w/2.)/w*2.],
        name = 'Shear Center',
        hovertext = r'Shear Center <br> x_sc = {0:.2f} <br> y_sc = {1:.2f} <br> w.r.t. ext. of bottom flange & web'.format(y_sc,x_sc),
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'blue', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the location of the center of gravity
    trace1_cg = go.Scatter(
        x = [origin_x + x_cdg],
        y = [origin_y + y_cdg],
        hovertext = r'Center of Gravity <br> x_cog = {0:.2f} <br> y_cog = {1:.2f} <br> w.r.t. ext. of bottom flange & web'.format(x_cdg,y_cdg),
        name = 'Center of Gravity',
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'purple', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the neutral axis
    phi = np.arctan((2*I_xy)/max([1e-05,(I_x - I_y)]))/2.

    trace1_na = go.Scatter(
        x = np.linspace(origin_x + x_cdg - max([h_l,h_r])/4.*np.cos(phi),origin_x + x_cdg + max([h_l,h_r])/4.*np.cos(phi),21),
        y = np.linspace(origin_y + y_cdg - max([h_l,h_r])/4.*np.sin(phi), origin_y + y_cdg + max([h_l,h_r])/4.*np.sin(phi), 21),
        name = 'Neutral Axis',
        mode='lines',
        hoveron = 'points',
        hoverinfo = 'name',
        line = dict(
            color = ('purple'),
            dash = 'dashdot'))

    #Determination of the shear flow fields in the upper and lower flanges of the I-section    
    x_par_flow_f = np.linspace(0.,w,21)
    y_par_flow_f = q23_y(x_par_flow_f[::-1])


    y_par_coord_f = [origin_x - np.abs(y_par_flow_f[x])*max([h_l,h_r])/5./max([q43_y(h_r),q12_y(h_l),max(q23_y(np.linspace(0.,w,21)))]) for x in range(len(y_par_flow_f))]
    y_par_flow_f = y_par_flow_f.tolist()
    y_par_flow_f.append(y_par_flow_f[-1])
    y_par_flow_f.insert(0,y_par_flow_f[0])



    # #Points through which the parabola-shaped shear flow in the web goes
    x_par_f = np.linspace(origin_y + t_w_l/2. ,origin_y + w - t_w_r/2.,21)

    # x_par_coord_l = x_par_coord_l.tolist()
    y_par_coord_f.append(origin_x)
    y_par_coord_f.insert(0,origin_x)


    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    x_par_f = x_par_f.tolist()
    x_par_f.append(origin_y + w - t_w_r/2.)
    x_par_f.insert(0,origin_y + t_w_l/2.)

    #Entity that plots the upper flange shear flow
    trace2_l = go.Scatter(
        x= y_par_coord_f,
        y= x_par_f,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in y_par_flow_f],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])


    #Computation of parabola constats for plotting the shear flow on web
    #Points through which the parabola-shaped shear flow in the web goes in order to define the appropriate tags

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow_l = np.linspace(0.,h_l,21)
    x_par_flow_l = q12_y(y_par_flow_l)

    x_par_coord_l = [origin_y - (x_par_flow_l[len(x_par_flow_l) - x -1])*max([h_l,h_r])/5./max([q43_y(h_r),q12_y(h_l),max(q23_y(np.linspace(0.,w,21)))]) for x in range(len(x_par_flow_l))]
    x_par_flow_l = x_par_flow_l.tolist()
    x_par_flow_l.append(x_par_flow_l[-1])
    x_par_flow_l.insert(0,x_par_flow_l[0])

    #Points through which the parabola-shaped shear flow in the web goes
    y_par_l = np.linspace(origin_x + t_f_u/2.,origin_x+h_l,21)

    #Points through which the parabola-shaped shear flow in the web goes
    x_par_coord_l.append(origin_y)
    x_par_coord_l.insert(0,origin_y)

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_l = y_par_l.tolist()
    y_par_l.append(origin_x + h_l)
    y_par_l.insert(0,origin_x + t_f_u/2.)

    #Entity that plots the web shear flow
    trace4_l = go.Scatter(
        x= y_par_l,
        y= x_par_coord_l,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in x_par_flow_l[::-1]],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow_r = np.linspace(0 ,h_r,21)
    x_par_flow_r = q43_y(y_par_flow_r)

    x_par_coord_r = [origin_y + w + (x_par_flow_r[len(x_par_flow_r) - x - 1])*max([h_l,h_r])/5./max([q43_y(h_r),q12_y(h_l),max(q23_y(np.linspace(0.,w,21)))]) for x in range(len(x_par_flow_r))]

    x_par_flow_r = x_par_flow_r.tolist()
    x_par_flow_r.append(x_par_flow_r[-1])
    x_par_flow_r.insert(0,0.)

    # #Points through which the parabola-shaped shear flow in the web goes
    y_par_r = np.linspace(origin_x + t_f_u/2.,origin_x+ h_r,21)

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_r = y_par_r.tolist()
    y_par_r.append(origin_x+h_r)
    y_par_r.insert(0,origin_x +t_f_u/2.)

    x_par_coord_r.append(origin_y + w)
    x_par_coord_r.insert(0,origin_y +  w)

    #Entity that plots the web shear flow
    trace4_r = go.Scatter(
        x= y_par_r,
        y= x_par_coord_r,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in x_par_flow_r[::-1]],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Set the origin of coordinates for the right-most figure (idealized)
    if max([h_r,h_l]) > w:
        origin_xd = (w/max([h_r,h_l])*(w - 0.75*max([h_r,h_l])) + 0.75*max([h_r,h_l]))
    else:
        origin_xd = 1.*w
    origin_yd = 0.

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = np.linspace(origin_yd  + t_w_l/2. ,origin_yd + w- t_w_r/2.,num+1).tolist()
    x_w_l =  np.linspace(origin_yd + t_w_l/2., origin_yd + t_w_l/2.,num+1).tolist()
    x_w_r =  np.linspace(origin_yd + w - t_w_r/2., origin_yd + w - t_w_r/2.,num+1).tolist()


    y_fu = np.linspace(origin_xd + t_f_u/2., origin_xd + t_f_u/2.,num+1).tolist()
    y_w_l =  np.linspace(origin_xd + t_f_u/2.,origin_xd + h_l,num+1).tolist()
    y_w_r =  np.linspace(origin_xd + t_f_u/2.,origin_xd + h_r,num+1).tolist()

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    x_w_l =list(itertools.chain.from_iterable(itertools.repeat(x_w_l[x], 2) for x in range(len(x_w_l))))
    x_w_r =list(itertools.chain.from_iterable(itertools.repeat(x_w_r[x], 2) for x in range(len(x_w_r))))

    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))
    y_w_l =list(itertools.chain.from_iterable(itertools.repeat(y_w_l[x], 2) for x in range(len(y_w_l))))
    y_w_r =list(itertools.chain.from_iterable(itertools.repeat(y_w_r[x], 2) for x in range(len(y_w_r))))
    
    x_fu_int = np.linspace(0. ,max([1,w]),num+1).tolist()
    x_fu_int = x_fu_int[::-1]
    y_w_l_int =  np.linspace(0.,max([1,h_l]),num+1).tolist()
    y_w_r_int =  np.linspace(0.,max([1,h_r]),num+1).tolist()

    x_fu_int = list(itertools.chain.from_iterable(itertools.repeat(x_fu_int[x], 2) for x in range(len(x_fu_int))))
    y_w_l_int =list(itertools.chain.from_iterable(itertools.repeat(y_w_l_int[x], 2) for x in range(len(y_w_l_int))))
    y_w_r_int =list(itertools.chain.from_iterable(itertools.repeat(y_w_r_int[x], 2) for x in range(len(y_w_r_int))))


    #Append constants at beginning and end of the list to offset shape as much as the shear flow value at the intersection with the following element
    y_flange_u = q23_y
    x_web_l = q12_y
    x_web_r = q43_y

    #Define discretized shear flow by averaging calculated flow on non-idealized cross section in flanges and web
    q_fu_s =  [integrate.quad(y_flange_u,x_fu_int[2*x],x_fu_int[2*x+2])[0]/(x_fu_int[2*x + 2] - x_fu_int[2*x]) for x in range(int(len(x_fu_int)/2 -1))]
    q_w_l_s =  [integrate.quad(x_web_l,y_w_l_int[2*x],y_w_l_int[2*x+2])[0]/(y_w_l_int[2*x + 2] - y_w_l_int[2*x]) for x in range(int(len(y_w_l_int)/2 -1))]
    q_w_r_s =  [integrate.quad(x_web_r,y_w_r_int[2*x],y_w_r_int[2*x+2])[0]/(y_w_r_int[2*x + 2] - y_w_r_int[2*x]) for x in range(int(len(y_w_r_int)/2 -1))]

    q_fu = [0.] + [q_fu_s[int(x/2)] for x in range(int(len(x_fu)-2))] + [0.]
    q_w_l = [0.] + [q_w_l_s[int(x/2)] for x in range(int(len(x_w_l)-2))] + [0.]
    q_w_r = [0.] + [q_w_r_s[int(x/2)] for x in range(int(len(x_w_r)-2))] + [0.]

    q_fu_pl = [np.abs(x) + origin_xd - max([h_l,h_r])/5.*np.abs(x)/max([q43_y(h_r),q12_y(h_l),max(q23_y(np.linspace(0.,w,21)))]) for x in q_fu]
    q_w_l_pl = [origin_yd - max([h_l,h_r])/5.*x/max([q43_y(h_r),q12_y(h_l),max(q23_y(np.linspace(0.,w,21)))]) for x in q_w_l[::-1]]
    q_w_r_pl = [origin_yd +w+ max([h_l,h_r])/5.*x/max([q43_y(h_r),q12_y(h_l),max(q23_y(np.linspace(0.,w,21)))]) for x in q_w_r[::-1]]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]
    q_w_l[0] = q_w_l_s[0]
    q_w_l[-1] = q_w_l_s[-1]
    q_w_r[0] = q_w_r_s[0]
    q_w_r[-1] = q_w_r_s[-1]

    #Calculate the corresponding boom areas based on the ratios of shear flow between two adjacent booms
    A_w_r = np.zeros(num+1)
    A_w_l = np.zeros(num+1)
    A_fu = np.zeros(num+1)
    
    y_w_l_A =  np.linspace(origin_xd + t_f_u/2.,origin_xd + h_l,num+1).tolist()
    y_w_r_A =  np.linspace(origin_xd + t_f_u/2.,origin_xd + h_r,num+1).tolist()
    x_fu_A = np.linspace(origin_yd  + t_w_l/2. ,origin_yd + w- t_w_r/2.,num+1).tolist()

    for x in range(0,num):
        A_w_r[num - x] += t_w_r*h_r/(num -1)/6.*(2 + q43_y(y_w_r_A[x+1] - origin_xd)/q43_y(y_w_r_A[x] - origin_xd))
        A_w_r[num - x -1] += t_w_r*h_r/(num -1)/6.*(2 + q43_y(y_w_r_A[x] - origin_xd)/q43_y(y_w_r_A[x+1] - origin_xd))

        A_w_l[num - x] += t_w_l*h_l/(num -1)/6.*(2 + q12_y(y_w_l_A[x+1] - origin_xd)/q12_y(y_w_l_A[x] - origin_xd))
        A_w_l[num - x -1] += t_w_l*h_l/(num -1)/6.*(2 + q12_y(y_w_l_A[x] - origin_xd)/q12_y(y_w_l_A[x+1] - origin_xd))

   
    for x in range(0, num):
        A_fu[x] += t_f_u*w/(num+1)/6.*(2 + q23_y(x_fu_A[x+1] - origin_yd)/q23_y(x_fu_A[x] - origin_yd))
        A_fu[x+1] += t_f_u*w/(num+1)/6.*(2 + q23_y(x_fu_A[x] - origin_yd)/q23_y(x_fu_A[x+1] - origin_yd))
    
    
    A_w_r[0] = A_w_r[0] + A_fu[-1]
    A_w_l[0] = A_w_l[0] + A_fu[0]
    A_fu[-1] = A_w_r[0]
    A_fu[0] = A_w_l[0]

    A_w_r = A_w_r.tolist()
    A_w_l = A_w_l.tolist()
    A_fu = A_fu.tolist()

    res_A = A_w_l + A_fu + A_w_r
    res_A = list(itertools.chain.from_iterable(itertools.repeat(res_A[x], 2) for x in range(len(res_A))))


    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x =  y_w_l + y_fu + y_w_r,
        y =  x_w_l + x_fu + x_w_r,
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['Ab = {0:.0f} mm2'.format(k) for k in res_A],
        line=dict(color='grey'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the upper flange
    trace2b = go.Scatter(
        x= q_fu_pl,
        y=  x_fu,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_fu],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the web
    trace4b_l = go.Scatter(
        x= y_w_l,
        y= q_w_l_pl,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_w_l[::-1]],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    trace4b_r = go.Scatter(
        x= y_w_r,
        y= q_w_r_pl,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_w_r[::-1]],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Add all the plotting entities to a data list that maps the plotting
    data_temp = [trace1, trace1_sc, trace1_cg, trace1_na, trace2_l, trace4_l, trace4_r, trace1b, trace2b, trace4b_l, trace4b_r]

    data = data + data_temp
    c = c +1

    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = [dict(x = origin_x + 8*h_r/10., y = origin_y + w - t_w_r/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=36, ay=0),
                       dict(x = origin_x + 4.5*h_r/10., y = origin_y + w - t_w_r/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=54, ay=0),
                       dict(x = origin_x + t_f_u, y = origin_y + w - t_w_r/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=66, ay=0),
                       
                       dict(x = origin_x + 0.98*h_l, y = origin_y + t_w_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-36, ay=0),
                       dict(x = origin_x + 0.7*h_l, y = origin_y + t_w_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-54, ay=0),
                       dict(x = origin_x + 0.36*h_l, y = origin_y + t_w_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-66, ay=0),

                       dict(x = origin_x + t_f_u/2., y = origin_y + t_w_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x + t_f_u/2., y = origin_y + 0.85*w, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x + t_f_u/2., y = origin_y + 0.2*w, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x + t_f_u/2., y = origin_y + 0.65*w, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x + t_f_u/2., y = origin_y + 0.37*w, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),

                       dict(x = origin_x - y_sc + t_f_u/2. , y =  origin_y + x_sc - t_w_r/2.*(x_sc - w/2.)/w*2. + 0.05*h_l, xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 14 ),showarrow = True, arrowhead=1, arrowsize=1.5, arrowwidth=1, arrowcolor='black', ax=0., ay= -40, xanchor = "center" ) ],
      #Define axes visibility properties
        xaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            # range=[-1.5*max([w,w_l]), 1.5*max([w,w_l])],
            showticklabels=False
        ),
        yaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            # range=[-0.2*h_l, 1.2*h_l],
            scaleanchor="x", 
            scaleratio=1,
            showticklabels=False
        )
    )

    return (data, layout)

def gen_data_s_beam(V_y,V_x,t_w_l,t_w_r,t_f_u,t_f_l,h,w_u,num):

    #Calculate geometrical properties of section
    A_tot = w_u*t_f_u + t_w_l*(h - t_f_u - t_f_l) + t_w_r*(h - t_f_u - t_f_l) + w_u*t_f_l #mm2
    y_cdg = (w_u*t_f_u*(h - t_f_u/2.) + t_w_l*(h - t_f_u - t_f_l)*h/2. + t_w_r*(h - t_f_u - t_f_l)*h/2. + w_u*t_f_l**2/2)/A_tot #mm w.r.t. bottom line
    x_cdg = (w_u*t_f_u*(w_u/2.) + t_w_l*(h - t_f_u - t_f_l)*t_w_l/2. + t_w_r*(h - t_f_u - t_f_l)*(w_u - t_w_r/2.) + w_u*t_f_l*(w_u/2.))/A_tot #mm w.r.t. ext left web
    I_x = 1./12.*w_u*t_f_u**3 + 1./12.*w_u*t_f_l**3 + w_u*t_f_u*(h - y_cdg - t_f_u/2.)**2 + w_u*t_f_l*(y_cdg - t_f_l/2.)**2 + 1./12.*t_w_l*(h - t_f_u - t_f_l)**3 +(h - t_f_u - t_f_l)*t_w_l*(h/2. - y_cdg)**2 + 1./12.*t_w_r*(h - t_f_u - t_f_l)**3 +(h - t_f_u - t_f_l)*t_w_r*(h/2. - y_cdg)**2  #mm4
    I_y = 1./12.*t_f_u*w_u**3 + 1./12.*t_f_l*w_u**3 + 1./12.*(h-t_f_u-t_f_l)*t_w_l**3 + (h-t_f_u-t_f_l)*t_w_l*(w_u/2. - t_w_l/2.)**2 + 1./12.*(h-t_f_u-t_f_l)*t_w_r**3 + (h-t_f_u-t_f_l)*t_w_r*(w_u/2. - t_w_r/2.)**2


    #Compute maximum shear flows on section
    q_flange_u = V_y*(t_f_u*w_u/2.*(h - y_cdg - t_f_u/2.))/I_x #N/mm
    q_flange_l = V_y*(t_f_l*w_u/2.*(y_cdg - t_f_l/2.))/I_x #N/mm

    Q_w_u = t_f_u*w_u*(h - y_cdg - t_f_u/2.) + (h - y_cdg - t_f_u)**2/2*t_w_l #mm3
    Q_w_u = t_f_l*w_u*(y_cdg - t_f_l/2.) + (y_cdg - t_f_l)**2/2*t_w_l #mm3

    q_web = V_y*Q_w_u/I_x #N/mm @ cog

    q_web_fun = lambda y: V_y*(t_f_u*w_u/2.*(h - y_cdg - t_f_u/2.) + (h - y - t_f_u)*t_w_l*(((h - t_f_u) + y)/2. - y_cdg))/I_x

    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    if max([w_u,w_u]) > h:
        origin_x = -(h/max([w_u,w_u])*(h - 0.75*max([w_u,w_u])) + 0.75*max([w_u,w_u]))
    else:
        origin_x = -1.*h
    origin_y = 0.

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1_o = go.Scatter(
        x=[origin_x - w_u/2., origin_x + w_u/2. ,origin_x+w_u/2., origin_x - w_u/2.,origin_x - w_u/2.],
        y=[origin_y,origin_y,origin_y+h,origin_y+h,origin_y],
        mode='lines',
        name=r'Rectangular Beam <br> V = {0:.2f} N <br> h = {1:.0f} mm <br> w = {2:.0f} mm <br> t1 = {3:.0f} mm <br> t2 = {4:.0f} mm <br> t3 = {5:.0f} mm <br> t4 = {6:.0f} mm <br> A = {7:.0f} mm2 <br> I_x = {8:.2f} mm4 <br> I_y = {9:.0f} mm4'.format(V_y,h,w_u,t_f_u,t_f_l,t_w_l,t_w_r,A_tot,I_x,I_y),
        hoverinfo='text',
        hoveron='fills',
        fill='tonexty', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    trace1_i = go.Scatter(
        x=[origin_x - w_u/2. + t_w_l, origin_x + w_u/2. - t_w_r ,origin_x+w_u/2. - t_w_r, origin_x - w_u/2. + t_w_l,origin_x - w_u/2. + t_w_l],
        y=[origin_y + t_f_l,origin_y + t_f_l ,origin_y+h - t_f_u,origin_y+h - t_f_u,origin_y + t_f_l],
        mode='lines',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    #Entity that plots the location of the shear load application (shear center)
    trace1_sc = go.Scatter(
        x = [origin_x],
        y = [y_cdg],
        name = 'Shear Center',
        hovertext = 'Center of Gravity & <br> Shear Center <br> x_cog = {0:.2f} <br> y_cog = {1:.2f} <br> w.r.t. ext. of bottom flange & left web'.format(x_cdg,y_cdg),
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'blue', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the neutral axis
    trace1_na = go.Scatter(
        x = np.linspace(origin_x - w_u/4.,origin_x+ w_u/4.,21),
        y = np.linspace(origin_y + y_cdg, origin_y + y_cdg, 21),
        name = 'Neutral Axis',
        mode='lines',
        hoveron = 'points',
        hoverinfo = 'name',
        line = dict(
            color = ('purple'),
            dash = 'dashdot'))

    #Determination of the shear flow fields in the upper and lower flanges of the I-section
    y_flange_u = np.linspace(q_flange_u,0.,11).tolist() + np.linspace(0.,q_flange_u,11).tolist()
    y_flange_l = np.linspace(q_flange_l,0.,11).tolist() + np.linspace(0.,q_flange_l,11).tolist()

    #Entity that plots the upper flange shear flow
    trace2 = go.Scatter(
        x= [origin_x - w_u/2. + t_w_l/2.] + np.linspace(origin_x - w_u/2 + t_w_l/2.,origin_x,11).tolist() + np.linspace(origin_x,origin_x + w_u/2 - t_w_r/2. ,11).tolist() + [origin_x + w_u/2. - t_w_r/2.],
        y= [origin_y + h] + np.linspace(origin_y + h + h/5.*q_flange_u/max([q_flange_u,q_flange_l]),origin_y + h,11).tolist() + np.linspace(origin_y + h,origin_y + h + h/5.*q_flange_u/max([q_flange_u,q_flange_l]),11).tolist() + [origin_y + h],
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
        x= [origin_x - w_u/2. + t_w_l/2.] + np.linspace(origin_x - w_u/2. + t_w_l/2.,origin_x,11).tolist() + np.linspace(origin_x,origin_x + w_u/2. - t_w_r/2.,11).tolist() + [origin_x + w_u/2. - t_w_r/2.],
        y= [origin_y] + np.linspace(origin_y - h/5.*q_flange_u/max([q_flange_u,q_flange_l]),origin_y,11).tolist() + np.linspace(origin_y ,origin_y - h/5.*q_flange_u/max([q_flange_u,q_flange_l]),11).tolist() + [origin_y],
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
    #Points through which the parabola-shaped shear flow in the web goes in order to define the appropriate tags
    e1_flow = [q_flange_l,t_f_l/2.]
    e2_flow = [q_flange_u,h - t_f_u/2.]
    vert_flow = [q_web, y_cdg]

    #Solve for the parabola constants
    a_flow, b_flow, c_flow = opt.fsolve(eq_parab,[1/180.,-5./18.,713./36.],args=(e1_flow,e2_flow,vert_flow))

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow = np.linspace(origin_y+t_f_l,origin_y+h-t_f_u,21)
    x_par_flow = q_web_fun(y_par_flow)

    x_par_coord_l = [origin_x - w_u/2. - x_par_flow[x]*h/5./max([q_flange_u,q_flange_l]) for x in range(len(x_par_flow))]
    x_par_coord_r = [origin_x + w_u/2. + x_par_flow[x]*h/5./max([q_flange_u,q_flange_l]) for x in range(len(x_par_flow))]
    x_par_flow = x_par_flow.tolist()
    x_par_flow.append(q_flange_u)
    x_par_flow.insert(0,q_flange_l)

    #Points through which the parabola-shaped shear flow in the web goes
    e1_coord = [origin_x - w_u/2. - h/5.*q_flange_l/max([q_flange_u,q_flange_l]),t_f_l/2.]
    e2_coord = [origin_x - w_u/2.- h/5.*q_flange_u/max([q_flange_u,q_flange_l]),h - t_f_u/2.]
    vert_coord = [origin_x -w_u/2.- h/5.*q_web/max([q_flange_u,q_flange_l]),y_cdg]

    #Solve for the parabola constants
    a_coord, b_coord, c_coord = opt.fsolve(eq_parab,[1/180.,-5./18.,713./36.],args=(e1_coord,e2_coord,vert_coord))

    y_par = np.linspace(origin_y+t_f_l/2.,origin_y+h-t_f_u/2.,21)

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par = y_par.tolist()
    y_par.append(origin_y+h-t_f_u/2.)
    y_par.insert(0,origin_y+t_f_l/2.)

    # x_par_coord = x_par_coord.tolist()
    x_par_coord_l.append(origin_x - w_u/2.)
    x_par_coord_l.insert(0,origin_x - w_u/2.)

    x_par_coord_r.append(origin_x + w_u/2.)
    x_par_coord_r.insert(0,origin_x + w_u/2.)

    #Entity that plots the web shear flow
    trace4_l = go.Scatter(
        x=x_par_coord_l,
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

    trace4_r = go.Scatter(
        x=x_par_coord_r,
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
    if max([w_u,w_u]) > h:
        origin_xd = (h/max([w_u,w_u])*(h - 0.75*max([w_u,w_u])) + 0.75*max([w_u,w_u]))
    else:
        origin_xd = 1.*h
    origin_yd = 0.

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = np.linspace(origin_xd - w_u/2. ,origin_xd+w_u/2.,num).tolist()
    x_fl = np.linspace(origin_xd - w_u/2. ,origin_xd+w_u/2.,num).tolist()
    x_w_l = np.linspace(origin_xd - w_u/2., origin_xd - w_u/2.,num+1).tolist()
    x_w_r = np.linspace(origin_xd + w_u/2., origin_xd + w_u/2.,num+1).tolist()


    y_fu = np.linspace(origin_yd + h - t_f_u/2., origin_yd + h - t_f_u/2.,num).tolist()
    y_fl = np.linspace(origin_yd + t_f_l/2., origin_yd + t_f_l/2.,num).tolist()
    y_w_l =  np.linspace(origin_yd + t_f_l/2.,origin_yd + h - t_f_u/2.,num+1).tolist()
    y_w_r =  np.linspace(origin_yd + t_f_l/2.,origin_yd + h - t_f_u/2.,num+1).tolist()

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    x_fl = list(itertools.chain.from_iterable(itertools.repeat(x_fl[x], 2) for x in range(len(x_fl))))
    x_w_l =list(itertools.chain.from_iterable(itertools.repeat(x_w_l[x], 2) for x in range(len(x_w_l))))
    x_w_r =list(itertools.chain.from_iterable(itertools.repeat(x_w_r[x], 2) for x in range(len(x_w_r))))

    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))
    y_fl = list(itertools.chain.from_iterable(itertools.repeat(y_fl[x], 2) for x in range(len(y_fl))))
    y_w_l =list(itertools.chain.from_iterable(itertools.repeat(y_w_l[x], 2) for x in range(len(y_w_l))))
    y_w_r =list(itertools.chain.from_iterable(itertools.repeat(y_w_r[x], 2) for x in range(len(y_w_r))))

    #Append constants at beginning and end of the list to offset shape as much as the shear flow value at the intersection with the following element
    y_flange_u = np.linspace(0.,q_flange_u,int(num/2 + 1)).tolist() + np.linspace(q_flange_u,0.,int(num/2 + 1)).tolist()[1:]
    y_flange_l = np.linspace(0.,q_flange_l,int(num/2 + 1)).tolist() + np.linspace(q_flange_l,0.,int(num/2 + 1)).tolist()[1:]

    x_web_l = q_web_fun
    x_web_r = q_web_fun


    #Define discretized shear flow by averaging calculated flow on non-idealized cross section in flanges and web
    q_fu_s = [(y_flange_u[x+1] + y_flange_u[x])/2 for x in range(len(y_flange_u) - 1)]
    q_fl_s = [(y_flange_l[x+1] + y_flange_l[x])/2 for x in range(len(y_flange_l) - 1)]
    q_wl_s =  [integrate.quad(x_web_l,y_w_l[2*x],y_w_l[2*x+2])[0]/(y_w_l[2*x + 2] - y_w_l[2*x]) for x in range(int(len(y_w_l)/2 -1))]
    q_wr_s =  [integrate.quad(x_web_r,y_w_r[2*x],y_w_r[2*x+2])[0]/(y_w_r[2*x + 2] - y_w_r[2*x]) for x in range(int(len(y_w_r)/2 -1))]

    q_w_l_total = integrate.quad(x_web_l,y_w_l[0],y_w_l[-1])[0]
    q_w_r_total = integrate.quad(x_web_r,y_w_r[0],y_w_r[-1])[0]

    q_fu = [0.] + [q_fu_s[int(x/2)] for x in range(int(len(x_fu)-2))] + [0.]
    q_fl = [0.] + [q_fl_s[int(x/2)] for x in range(int(len(x_fl)-2))] + [0.]
    q_wl = [0.] + [q_wl_s[int(x/2)] for x in range(int(len(x_w_l)-2))] + [0.]
    q_wr = [0.] + [q_wr_s[int(x/2)] for x in range(int(len(x_w_r)-2))] + [0.]

    q_fu_pl = [x + origin_y + h + h/5.*x/max([q_flange_u,q_flange_l]) for x in q_fu]
    q_fl_pl = [x + origin_y  - h/5.*x/max([q_flange_u,q_flange_l]) for x in q_fl]
    q_wl_pl = [origin_xd -w_u/2. - h/5.*x/max([q_flange_u,q_flange_l]) for x in q_wl]
    q_wr_pl = [origin_xd +w_u/2. + h/5.*x/max([q_flange_u,q_flange_l]) for x in q_wr]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]
    q_fl[0] = q_fl_s[0]
    q_fl[-1] = q_fl_s[-1]
    q_wl[0] = q_wl_s[0]
    q_wl[-1] = q_wl_s[-1]
    q_wr[0] = q_wr_s[0]
    q_wr[-1] = q_wr_s[-1]

    #Calculate the corresponding boom areas based on the ratios of shear flow between two adjacent booms
    A_fu = np.zeros(num)
    A_fl = np.zeros(num)
    A_wl = np.zeros(num+1)
    A_wr = np.zeros(num+1)


    y_w_A =  np.linspace(origin_yd,origin_yd + h - t_f_u/2.,num+1).tolist()
    x_fu_A = np.linspace(origin_xd ,origin_xd+w_u/2.,num).tolist()
    x_fl_A = np.linspace(origin_xd ,origin_xd+w_u/2.,num).tolist()

    q_fupper = lambda s: q_flange_u*( 1 -1/w_u*2*s)

    for x in range(0,math.ceil(num/2) - 1):
        A_fu[x] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x] - origin_xd)]))
        A_fu[x +1] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x +1] - origin_xd)]))
        
        A_fu[num - x -1] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x] - origin_xd)]))
        A_fu[num - x -2] += t_f_u*w_u/(num -1)/6.*(2 + q_fupper(x_fu_A[x] - origin_xd)/max([1e-05,q_fupper(x_fu_A[x +1] - origin_xd)]))
        
        A_fl[x] += t_f_l*w_u/(num -1)/6.*(2 + q_fupper(x_fl_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x] - origin_xd)]))
        A_fl[x +1] += t_f_l*w_u/(num -1)/6.*(2 + q_fupper(x_fl_A[x] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x +1] - origin_xd)]))
        
        A_fl[num - x -1] += t_f_l*w_u/(num -1)/6.*(2 + q_fupper(x_fl_A[x+1] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x] - origin_xd)]))
        A_fl[num - x -2] += t_f_l*w_u/(num -1)/6.*(2 + q_fupper(x_fl_A[x] - origin_xd)/max([1e-05,q_fupper(x_fl_A[x +1] - origin_xd)]))


    for x in range(0, num):
        A_wl[x] += t_w_l*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x+1] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x] - origin_yd)]))
        A_wl[x+1] += t_w_l*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x+1] - origin_yd)]))
        
        A_wr[x] += t_w_r*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x+1] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x] - origin_yd)]))
        A_wr[x+1] += t_w_r*h/(num+1)/6.*(2 + q_web_fun(y_w_A[x] - origin_yd)/max([1e-05,q_web_fun(y_w_A[x+1] - origin_yd)]))


    A_fu[0] = A_fu[0] + A_wl[-1]
    A_fl[0] = A_fl[0] + A_wl[0]
    A_wl[-1] = A_fu[0]
    A_wl[0] = A_fl[0]

    A_fu[-1] = A_fu[-1] + A_wr[-1]
    A_fl[-1] = A_fl[-1] + A_wr[0]
    A_wr[-1] = A_fu[-1]
    A_wr[0] = A_fl[-1]

    A_fu = A_fu.tolist()
    A_fl = A_fl.tolist()
    A_wl = A_wl.tolist()
    A_wr = A_wr.tolist()

    res_A = A_fl + A_wl + A_fu + A_wr
    res_A = list(itertools.chain.from_iterable(itertools.repeat(res_A[x], 2) for x in range(len(res_A))))


    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x = x_fl + x_w_l + x_fu + x_w_r,
        y = y_fl + y_w_l + y_fu + y_w_r,
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['Ab = {0:.0f} mm2'.format(k) for k in res_A],
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
    trace4b_l = go.Scatter(
        x=q_wl_pl,
        y=y_w_l,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_wl],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])
        #Entity that plots the discretized shear flow on the web
    
    trace4b_r = go.Scatter(
        x=q_wr_pl,
        y=y_w_r,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(k) for k in q_wr],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Add all the plotting entities to a data list that maps the plotting
    data_temp = [trace1_o, trace1_i, trace1_sc, trace1_na, trace2, trace3, trace4_l, trace4_r, trace1b, trace2b, trace3b, trace4b_l, trace4b_r]

    data = data + data_temp
    c = c +1

    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = [dict(x = origin_x + w_u/10.,y = origin_y + h - t_f_u/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-18,  ay=0),
                       dict(x = origin_x + w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-30, ay=0),
                       dict(x = origin_x + w_u/2. - t_w_l/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-45, ay=0),
                       dict(x = origin_x - w_u/10., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=18, ay=0),
                       dict(x = origin_x -w_u/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=30, ay=0),
                       dict(x = origin_x - w_u/2. + t_w_l/4., y = origin_y + h - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=45, ay=0),
                       
                       dict(x = origin_x - w_u/2. + 0.23*w_u,y = origin_y + t_f_l/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-45,  ay=0),
                       dict(x = origin_x - w_u/2. + 0.38*w_u, y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-30, ay=0),
                       dict(x = origin_x - w_u/2. + 0.49*w_u - t_w_l/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-18, ay=0),
                       dict(x = origin_x - w_u/2. +w_u - 0.23*w_u, y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=45, ay=0),
                       dict(x = origin_x - w_u/2. + w_u - 0.38*w_u, y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=30, ay=0),
                       dict(x = origin_x - w_u/2. + 0.49*w_u + t_w_l/4., y = origin_y + t_f_l/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=18, ay=0),

                       dict(x = origin_x - w_u/2. + t_w_l/2., y = origin_y + t_f_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x - w_u/2. + t_w_l/2., y = origin_y + 0.85*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x - w_u/2. + t_w_l/2., y = origin_y + 0.2*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x - w_u/2. + t_w_l/2., y = origin_y + 0.65*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x - w_u/2. + t_w_l/2., y = origin_y + 0.37*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),

                       dict(x = origin_x + w_u/2 - t_w_r/2., y = origin_y + t_f_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x + w_u/2 - t_w_r/2., y = origin_y + 0.85*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x + w_u/2 - t_w_r/2., y = origin_y + 0.2*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x + w_u/2 - t_w_r/2., y = origin_y + 0.65*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x + w_u/2 - t_w_r/2., y = origin_y + 0.37*h, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),
                       
                       dict(x = origin_x , y = y_cdg + 0.05*h, xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 14 ),showarrow = True, arrowhead=1, arrowsize=1.5, arrowwidth=1, arrowcolor='black', ax=0., ay= -40, xanchor = "center" ) ],
        #Define axes visibility properties
        xaxis=dict(
            autorange=False,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            range=[-1.5*max([w_u,w_u]), 1.5*max([w_u,w_u])],
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

def gen_data_u_beam(V_x,V_y,t_f_u,t_w_l,t_w_r,h_l,h_r,w,num):
    
    #Calculate geometrical properties of section
    A_tot = w*t_f_u + t_w_l*(h_l - t_f_u) + t_w_r*(h_r - t_f_u)
    y_cdg = (t_f_u*w*w/2. + t_w_r*h_r*w)/(t_f_u*w + t_w_r*h_r + t_w_l*h_l)
    x_cdg = (t_w_l*h_l*h_l/2. + t_w_r*h_r*h_r/2.)/(t_f_u*w + t_w_r*h_r + t_w_l*h_l)

    I_x = t_w_l*h_l*y_cdg**2 + (t_f_u*w**3/12) + t_f_u*w*(w/2. - y_cdg)**2 + t_w_r*h_r*(w - y_cdg)**2
    I_y = (t_w_l*h_l**3/12.) + t_w_l*h_l*(h_l/2. - x_cdg)**2 + t_f_u*w*x_cdg**2 + (t_w_r*h_r**3/12.) + t_w_r*h_r*(h_r/2. - x_cdg)**2
    I_xy = t_w_l*h_l*(h_l/2. - x_cdg)*y_cdg + t_f_u*w*(-w/2. + y_cdg)*(-x_cdg) + t_w_r*h_r*(h_r/2. - x_cdg)*(-w + y_cdg)
    
    #Compute maximum shear flows on section
    q43_y = lambda s: V_y*(I_xy/(I_x*I_y - I_xy**2)*t_w_r*((h_r - x_cdg)*s - s**2/2.) + I_y/(I_x*I_y - I_xy**2)*t_w_r*(w - y_cdg)*s)
    q43_x = lambda s: V_x*(I_x/(I_x*I_y - I_xy**2)*t_w_r*((h_r - x_cdg)*s - s**2/2.) + I_xy/(I_x*I_y - I_xy**2)*t_w_r*(w - y_cdg)*s)

    y_sc = integrate.quad(q43_y,0.,h_r)[0]*w/V_y
    x_sc = integrate.quad(q43_x,0.,h_r)[0]*w/V_x


    q12_y = lambda s: V_y*(-I_xy/(I_x*I_y - I_xy**2)*t_w_l*((h_l - x_cdg)*s - s**2/2.) + I_y/(I_x*I_y - I_xy**2)*t_w_l*(y_cdg)*s)
    q12_x = lambda s: V_x*(-I_x/(I_x*I_y - I_xy**2)*t_w_l*((h_l - x_cdg)*s - s**2/2.) + I_xy/(I_x*I_y - I_xy**2)*t_w_l*(y_cdg)*s)

    q23_x = lambda s: V_x*(I_x/(I_x*I_y - I_xy**2)*t_f_u*(-(x_cdg)*s) - I_xy/(I_x*I_y - I_xy**2)*t_f_u*(-(w - y_cdg)*s + s**2/2.)) + q43_x(h_r)


    q_web_fun_l = lambda y: V_l*(t_f_u*(w/2. + x_sc)*(h_l - y_cdg_l - t_f_u/2.) + (h_l - y - t_f_u)*t_w_l*(((h_l - t_f_u) + y)/2. - y_cdg_l))/I_x_l

    q_web_fun_r = lambda y: V_r*(t_f_u*(w/2. - x_sc)*(h_r - y_cdg_r - t_f_u/2.) + (h_r - y - t_f_u)*t_w_r*(((h_r - t_f_u) + y)/2. - y_cdg_r))/I_x_r

    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    if w > max([h_l,h_r]):
        origin_x = -(max([h_l,h_r])/w*(max([h_l,h_r]) - w) + w)
    else:
        origin_x = -1.*max([h_l,h_r])
    origin_y = 0.

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1 = go.Scatter(
        x=[origin_x -w/2. + t_w_l,origin_x -w/2. + t_w_l,origin_x+w/2. - t_w_r ,origin_x+w/2. - t_w_r ,origin_x + w/2.,origin_x + w/2., origin_x -w/2., origin_x - w/2., origin_x -w/2. + t_w_l],
        y=[origin_y + max([h_l,h_r]) - h_l,origin_y+max([h_l,h_r]) -t_f_u,origin_y+max([h_l,h_r]) -t_f_u,origin_y+max([h_l,h_r]) - h_r,origin_y+max([h_l,h_r]) - h_r,origin_y + max([h_l,h_r]),origin_y + max([h_l,h_r]), origin_y + max([h_l,h_r]) - h_l, origin_y + max([h_l,h_r]) - h_l],
        mode='lines',
        name=r'U Beam <br> V = {0:.2f} N <br> h1 = {1:.0f} mm <br> h2 = {2:.0f} mm <br> w = {3:.0f} mm <br> t2 = {4:.0f} mm <br> t3 = {5:.0f} mm <br> t4 = {6:.0f} mm <br> A = {7:.0f} mm2 <br> I_x = {8:.0f} mm4 <br> I_y = {9:.0f} mm4 <br> I_xy = {10:.0f} mm4'.format(V_y,h_l,h_r,w,t_f_u,t_w_l,t_w_r,A_tot,I_y,I_x,I_xy),
        hoverinfo='text',
        hoveron='fills',
        fill='toself', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    #Entity that plots the location of the shear load application (shear center)
    trace1_sc = go.Scatter(
        x = [origin_x - w/2. + x_sc - t_w_r/2.*(x_sc - w/2.)/w*2.],
        y = [origin_y + max([h_l,h_r]) + y_sc - t_f_u/2.],
        name = 'Shear Center',
        hovertext = r'Shear Center <br> x_sc = {0:.2f} <br> y_sc = {1:.2f} <br> w.r.t. ext. of upper flange & left web'.format(x_sc,y_sc),
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'blue', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the location of the center of gravity
    trace1_cg = go.Scatter(
        x = [origin_x - w/2. + y_cdg],
        y = [origin_y + max([h_l,h_r]) - x_cdg],
        hovertext = r'Center of Gravity <br> x_cog = {0:.2f} <br> y_cog = {1:.2f} <br> w.r.t. ext. of upper flange & left web'.format(y_cdg,x_cdg),
        name = 'Center of Gravity',
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'purple', size = 10),
        visible=[True if c is 0 else False][0])
   
    #Entity that plots the neutral axis
    phi = np.arctan(-(2*I_xy)/max([1e-05,(I_x - I_y)]))/2.

    trace1_na = go.Scatter(
        x = np.linspace(origin_x -w/2. + y_cdg - w/4.*np.cos(phi),origin_x -w/2. + y_cdg + w/4.*np.cos(phi),21),
        y = np.linspace(origin_y + max([h_l,h_r]) - x_cdg - w/4.*np.sin(phi), origin_y + max([h_l,h_r]) - x_cdg + w/4.*np.sin(phi), 21),
        name = 'Neutral Axis',
        mode='lines',
        hoveron = 'points',
        hoverinfo = 'name',
        line = dict(
            color = ('purple'),
            dash = 'dashdot'))


    #Determination of the shear flow fields in the upper and lower flanges of the I-section    
    x_par_flow_f = np.linspace(0.,w,21)
    y_par_flow_f = q23_x(x_par_flow_f)

    y_par_coord_f = [origin_y + max([h_l,h_r]) + np.abs(y_par_flow_f[x])*max([h_l,h_r])/5./max([q43_x(h_r),q12_x(h_l),max(q23_x(np.linspace(0.,w,21)))]) for x in range(len(y_par_flow_f))]
    y_par_flow_f = y_par_flow_f.tolist()
    y_par_flow_f.append(y_par_flow_f[-1])
    y_par_flow_f.insert(0,y_par_flow_f[0])

    # #Points through which the parabola-shaped shear flow in the web goes
    x_par_f = np.linspace(origin_x + w/2. - t_w_r ,origin_x - w/2. + t_w_l,21)

    # x_par_coord_l = x_par_coord_l.tolist()
    y_par_coord_f.append(origin_y + max([h_l,h_r]))
    y_par_coord_f.insert(0,origin_y + max([h_l,h_r]))


    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    x_par_f = x_par_f.tolist()
    x_par_f.append(origin_x - w/2. + t_w_l)
    x_par_f.insert(0,origin_x + w/2. - t_w_r)

    #Entity that plots the upper flange shear flow
    trace2_l = go.Scatter(
        x= x_par_f,
        y= y_par_coord_f,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in y_par_flow_f],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])


    #Computation of parabola constats for plotting the shear flow on web
    #Points through which the parabola-shaped shear flow in the web goes in order to define the appropriate tags

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow_l = np.linspace(0.,h_l,21)
    x_par_flow_l = q12_x(y_par_flow_l)

    x_par_coord_l = [origin_x - w/2. - np.abs(x_par_flow_l[x])*max([h_l,h_r])/5./max([q43_x(h_r),q12_x(h_l),max(q23_x(np.linspace(0.,w,21)))]) for x in range(len(x_par_flow_l))]
    x_par_flow_l = x_par_flow_l.tolist()
    x_par_flow_l.append(x_par_flow_l[-1])
    x_par_flow_l.insert(0,0.)

    # #Points through which the parabola-shaped shear flow in the web goes
    y_par_l = np.linspace(origin_y + max([h_l,h_r]) - h_l,origin_y+max([h_l,h_r])-t_f_u/2.,21)

    #Points through which the parabola-shaped shear flow in the web goes
    x_par_coord_l.append(origin_x - w/2.)
    x_par_coord_l.insert(0,origin_x -  w/2.)

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_l = y_par_l.tolist()
    y_par_l.append(origin_y+max([h_l,h_r])-t_f_u/2.)
    y_par_l.insert(0,origin_y + max([h_l,h_r]) - h_l)

    #Entity that plots the web shear flow
    trace4_l = go.Scatter(
        x=x_par_coord_l,
        y=y_par_l,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in x_par_flow_l],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_flow_r = np.linspace(0 ,h_r,21)
    x_par_flow_r = q43_x(y_par_flow_r)
    # print(x_par_flow_r)

    x_par_coord_r = [origin_x + w/2. + np.abs(x_par_flow_r[x])*max([h_l,h_r])/5./max([q43_x(h_r),q12_x(h_l),max(q23_x(np.linspace(0.,w,21)))]) for x in range(len(x_par_flow_r))]

    x_par_flow_r = x_par_flow_r.tolist()
    x_par_flow_r.append(x_par_flow_r[-1])
    x_par_flow_r.insert(0,0.)

    # #Points through which the parabola-shaped shear flow in the web goes
    y_par_r = np.linspace(origin_y + max([h_l,h_r]) - h_r,origin_y+max([h_l,h_r])-t_f_u/2.,21)

    #Append constants at beginning and end of the list to offset the parabola as much as the shear flow value at the intersection with the flanges
    y_par_r = y_par_r.tolist()
    y_par_r.append(origin_y+max([h_l,h_r])-t_f_u/2.)
    y_par_r.insert(0,origin_y +max([h_l,h_r]) - h_r)

    # x_par_coord = x_par_coord.tolist()
    x_par_coord_r.append(origin_x + w/2.)
    x_par_coord_r.insert(0,origin_x +  w/2.)

    #Entity that plots the web shear flow
    trace4_r = go.Scatter(
        x=x_par_coord_r,
        y=y_par_r,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in x_par_flow_r],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    #Set the origin of coordinates for the right-most figure (idealized)
    if w > max([h_l,h_r]):
        origin_xd = (max([h_l,h_r])/w*(max([h_l,h_r]) - w) + w)
    else:
        origin_xd = 1.*max([h_l,h_r])

    origin_yd = 0.

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = np.linspace(origin_xd - w/2. + t_w_l/2. ,origin_xd+w/2. - t_w_r/2.,num).tolist()
    x_w_l =  np.linspace(origin_xd - w/2. + t_w_l/2., origin_xd - w/2. + t_w_l/2.,num+1).tolist()
    x_w_r =  np.linspace(origin_xd + w/2. - t_w_r/2., origin_xd + w/2. - t_w_r/2.,num+1).tolist()


    y_fu = np.linspace(origin_yd + max([h_l,h_r]) - t_f_u/2., origin_yd + max([h_l,h_r]) - t_f_u/2.,num).tolist()
    y_w_l =  np.linspace(origin_yd + max([h_l,h_r]) - h_l,origin_yd + max([h_l,h_r]) - t_f_u/2.,num+1).tolist()
    y_w_r =  np.linspace(origin_yd + max([h_l,h_r]) - t_f_u/2.,origin_yd + max([h_l,h_r]) - h_r,num+1).tolist()

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    x_w_l =list(itertools.chain.from_iterable(itertools.repeat(x_w_l[x], 2) for x in range(len(x_w_l))))
    x_w_r =list(itertools.chain.from_iterable(itertools.repeat(x_w_r[x], 2) for x in range(len(x_w_r))))

    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))
    y_w_l =list(itertools.chain.from_iterable(itertools.repeat(y_w_l[x], 2) for x in range(len(y_w_l))))
    y_w_r =list(itertools.chain.from_iterable(itertools.repeat(y_w_r[x], 2) for x in range(len(y_w_r))))
    
    x_fu_int = np.linspace(0. ,max([1,w]),num).tolist()
    y_w_l_int =  np.linspace(0.,max([1,h_l]),num+1).tolist()
    y_w_r_int =  np.linspace(0.,max([1,h_r]),num+1).tolist()

    x_fu_int = list(itertools.chain.from_iterable(itertools.repeat(x_fu_int[x], 2) for x in range(len(x_fu_int))))
    y_w_l_int =list(itertools.chain.from_iterable(itertools.repeat(y_w_l_int[x], 2) for x in range(len(y_w_l_int))))
    y_w_r_int =list(itertools.chain.from_iterable(itertools.repeat(y_w_r_int[x], 2) for x in range(len(y_w_r_int))))


    #Append constants at beginning and end of the list to offset shape as much as the shear flow value at the intersection with the following element
    y_flange_u = q23_x
    x_web_l = q12_x
    x_web_r = q43_x

    #Define discretized shear flow by averaging calculated flow on non-idealized cross section in flanges and web
    q_fu_s =  [integrate.quad(y_flange_u,x_fu_int[2*x],x_fu_int[2*x+2])[0]/(x_fu_int[2*x + 2] - x_fu_int[2*x]) for x in range(int(len(x_fu_int)/2 -1))]
    q_w_l_s =  [integrate.quad(x_web_l,y_w_l_int[2*x],y_w_l_int[2*x+2])[0]/(y_w_l_int[2*x + 2] - y_w_l_int[2*x]) for x in range(int(len(y_w_l_int)/2 -1))]
    q_w_r_s =  [integrate.quad(x_web_r,y_w_r_int[2*x],y_w_r_int[2*x+2])[0]/(y_w_r_int[2*x + 2] - y_w_r_int[2*x]) for x in range(int(len(y_w_r_int)/2 -1))]



    q_fu = [0.] + [q_fu_s[int(x/2)] for x in range(int(len(x_fu)-2))] + [0.]
    q_w_l = [0.] + [np.abs(q_w_l_s[int(x/2)]) for x in range(int(len(x_w_l)-2))] + [0.]
    q_w_r = [0.] + [q_w_r_s[int(x/2)] for x in range(int(len(x_w_r)-2))] + [0.]



    q_fu_pl = [np.abs(x) + origin_y + max([h_l,h_r]) + max([h_l,h_r])/5.*np.abs(x)/max([q43_x(h_r),q12_x(h_l),max(q23_x(np.linspace(0.,w,21)))]) for x in q_fu[::-1]]
    q_w_l_pl = [origin_xd -w/2.- max([h_l,h_r])/5.*np.abs(x)/max([q43_x(h_r),q12_x(h_l),max(q23_x(np.linspace(0.,w,21)))]) for x in q_w_l]
    q_w_r_pl = [origin_xd +w/2.+ max([h_l,h_r])/5.*np.abs(x)/max([q43_x(h_r),q12_x(h_l),max(q23_x(np.linspace(0.,w,21)))]) for x in q_w_r[::-1]]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]
    q_w_l[0] = q_w_l_s[0]
    q_w_l[-1] = q_w_l_s[-1]
    q_w_r[0] = q_w_r_s[0]
    q_w_r[-1] = q_w_r_s[-1]

    #Calculate the corresponding boom areas based on the ratios of shear flow between two adjacent booms
    A_w_r = np.zeros(num+1)
    A_w_l = np.zeros(num+1)
    A_fu = np.zeros(num)
    
    y_w_l_A =  np.linspace(origin_yd,origin_yd  + h_l,num+1).tolist()
    y_w_r_A =  np.linspace(origin_yd,origin_yd  + h_r,num+1).tolist()
    x_fu_A = np.linspace(origin_xd ,origin_xd+w,num).tolist()

    for x in range(0,num):
        A_w_r[x] += t_w_r*h_r/(num -1)/6.*(2 + q43_x(y_w_r_A[x+1] - origin_yd)/max([1e-05,q43_x(y_w_r_A[x] - origin_yd)]))
        A_w_r[x +1] += t_w_r*h_r/(num -1)/6.*(2 + q43_x(y_w_r_A[x] - origin_yd)/max([1e-05,q43_x(y_w_r_A[x+1] - origin_yd)]))
        
        A_w_l[x] += t_w_l*h_l/(num -1)/6.*(2 + q12_x(y_w_l_A[x+1] - origin_yd)/min([-1e-05,q12_x(y_w_l_A[x] - origin_yd)]))
        A_w_l[x +1] += t_w_l*h_l/(num -1)/6.*(2 + q12_x(y_w_l_A[x] - origin_yd)/min([-1e-05,q12_x(y_w_l_A[x+1] - origin_yd)]))

    for x in range(0, num-1):
        A_fu[x] += t_f_u*w/(num+1)/6.*(2 + q23_x(x_fu_A[x+1] - origin_xd)/q23_x(x_fu_A[x] - origin_xd))
        A_fu[x+1] += t_f_u*w/(num+1)/6.*(2 + q23_x(x_fu_A[x] - origin_xd)/q23_x(x_fu_A[x+1] - origin_xd))
    
    # print(A_fu)
    A_w_r[-1] = A_w_r[-1] + A_fu[-1]
    A_w_l[-1] = A_w_l[-1] + A_fu[0]
    A_fu[-1] = A_w_r[-1]
    A_fu[0] = A_w_l[-1]

    A_w_r = A_w_r.tolist()
    A_w_l = A_w_l.tolist()
    A_fu = A_fu.tolist()

    res_A = A_w_l + A_fu + A_w_r
    res_A = list(itertools.chain.from_iterable(itertools.repeat(res_A[x], 2) for x in range(len(res_A))))


    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x =  x_w_l + x_fu + x_w_r,
        y =  y_w_l + y_fu + y_w_r[::-1],
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['Ab = {0:.0f} mm2'.format(k) for k in res_A],
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
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in q_fu],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the web
    trace4b_l = go.Scatter(
        x=q_w_l_pl,
        y=y_w_l,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in q_w_l],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    trace4b_r = go.Scatter(
        x=q_w_r_pl,
        y=y_w_r,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in q_w_r[::-1]],
        #hovertemplate = '%{hovertext:.2f}',
        fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    #Add all the plotting entities to a data list that maps the plotting
    data_temp = [trace1, trace1_cg, trace1_na, trace2_l, trace4_l, trace4_r, trace1_sc, trace1b, trace2b, trace4b_l, trace4b_r]

    data = data + data_temp
    c = c +1

    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = [dict(x = origin_x + w/10.,y = origin_y + max([h_l,h_r]) - t_f_u/2.,xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-18,  ay=0),
                       dict(x = origin_x + w/4., y = origin_y + max([h_l,h_r]) - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-30, ay=0),
                       dict(x = origin_x + w/2. - t_w_l/4., y = origin_y + max([h_l,h_r]) - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=-45, ay=0),
                       dict(x = origin_x - w/10., y = origin_y + max([h_l,h_r]) - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=18, ay=0),
                       dict(x = origin_x -w/4., y = origin_y + max([h_l,h_r]) - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=30, ay=0),
                       dict(x = origin_x - w/2. + t_w_l/4., y = origin_y + max([h_l,h_r]) - t_f_u/2., xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=45, ay=0),
                      
                       dict(x = origin_x - w/2. + t_w_l/2., y = origin_y+ max([h_l,h_r]) - h_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x - w/2. + t_w_l/2., y = origin_y+ max([h_l,h_r]) - h_l + 0.85*h_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x - w/2. + t_w_l/2., y = origin_y+ max([h_l,h_r]) - h_l + 0.2*h_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x - w/2. + t_w_l/2., y = origin_y+ max([h_l,h_r]) - h_l + 0.65*h_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x - w/2. + t_w_l/2., y = origin_y+ max([h_l,h_r]) - h_l + 0.37*h_l, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),

                       dict(x = origin_x + w/2. - t_w_r/2., y = origin_y + max([h_l,h_r]) - h_r, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x + w/2. - t_w_r/2., y = origin_y + max([h_l,h_r]) - h_r + 0.85*h_r, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -30 ),
                       dict(x = origin_x + w/2. - t_w_r/2., y = origin_y + max([h_l,h_r]) - h_r + 0.2*h_r, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x + w/2. - t_w_r/2., y = origin_y + max([h_l,h_r]) - h_r + 0.65*h_r, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -36 ),
                       dict(x = origin_x + w/2. - t_w_r/2., y = origin_y + max([h_l,h_r]) - h_r + 0.37*h_r, xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax=0., ay= -54 ),


                       dict(x = origin_x - w/2. + x_sc - t_w_r/2.*(x_sc - w/2.)/w*2., y = max([h_l,h_r]) + y_sc + 0.05*h_l, xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 14 ),showarrow = True, arrowhead=1, arrowsize=1.5, arrowwidth=1, arrowcolor='black', ax=0., ay= -40, xanchor = "center" ) ],
      #Define axes visibility properties
        xaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            # range=[-1.5*max([w,w_l]), 1.5*max([w,w_l])],
            showticklabels=False
        ),
        yaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            # range=[-0.2*h_l, 1.2*h_l],
            scaleanchor="x", 
            scaleratio=1,
            showticklabels=False
        )
    )

    return (data, layout)

def gen_data_o_beam(V_y,V_x,rad,t,num):
    

    #Calculate geometrical properties of section
    w = rad*2
    h_l = rad*2
    h_r = rad*2

    t_f_u = t
    t_w_l = t
    t_w_r = t


    A_tot = w*t_f_u + t_w_l*(h_l - t_f_u) + t_w_r*(h_r - t_f_u)
    y_cdg = (t_f_u*w*w/2. + t_w_r*h_r*w)/(t_f_u*w + t_w_r*h_r + t_w_l*h_l)
    x_cdg = (t_w_l*h_l*h_l/2. + t_w_r*h_r*h_r/2.)/(t_f_u*w + t_w_r*h_r + t_w_l*h_l)

    I_x = t_w_l*h_l*y_cdg**2 + (t_f_u*w**3/12) + t_f_u*w*(w/2. - y_cdg)**2 + t_w_r*h_r*(w - y_cdg)**2
    I_y = (t_w_l*h_l**3/12.) + t_w_l*h_l*(h_l/2. - x_cdg)**2 + t_f_u*w*x_cdg**2 + (t_w_r*h_r**3/12.) + t_w_r*h_r*(h_r/2. - x_cdg)**2
    I_xy = t_w_l*h_l*(h_l/2. - x_cdg)*y_cdg + t_f_u*w*(-w/2. + y_cdg)*(-x_cdg) + t_w_r*h_r*(h_r/2. - x_cdg)*(-w + y_cdg)
    

    #Compute maximum shear flows on section
    q43_y = lambda s: V_y*(I_xy/(I_x*I_y - I_xy**2)*t_w_r*((h_r - x_cdg)*s - s**2/2.) + I_y/(I_x*I_y - I_xy**2)*t_w_r*(w - y_cdg)*s)
    q43_x = lambda s: V_x*(I_x/(I_x*I_y - I_xy**2)*t_w_r*((h_r - x_cdg)*s - s**2/2.) + I_xy/(I_x*I_y - I_xy**2)*t_w_r*(w - y_cdg)*s)

    y_sc = integrate.quad(q43_y,0.,h_r)[0]*w/V_y
    x_sc = integrate.quad(q43_x,0.,h_r)[0]*w/V_x


    q12_y = lambda s: V_y*(-I_xy/(I_x*I_y - I_xy**2)*t_w_l*((h_l - x_cdg)*s - s**2/2.) + I_y/(I_x*I_y - I_xy**2)*t_w_l*(y_cdg)*s)
    q12_x = lambda s: V_x*(-I_x/(I_x*I_y - I_xy**2)*t_w_l*((h_l - x_cdg)*s - s**2/2.) + I_xy/(I_x*I_y - I_xy**2)*t_w_l*(y_cdg)*s)

    q23_x = lambda s: V_x*(I_x/(I_x*I_y - I_xy**2)*t_f_u*(-(x_cdg)*s) - I_xy/(I_x*I_y - I_xy**2)*t_f_u*(-(w - y_cdg)*s + s**2/2.)) + q43_x(h_r)

    q_web_fun_l = lambda y: V_l*(t_f_u*(w/2. + x_sc)*(h_l - y_cdg_l - t_f_u/2.) + (h_l - y - t_f_u)*t_w_l*(((h_l - t_f_u) + y)/2. - y_cdg_l))/I_x_l

    q_web_fun_r = lambda y: V_r*(t_f_u*(w/2. - x_sc)*(h_r - y_cdg_r - t_f_u/2.) + (h_r - y - t_f_u)*t_w_r*(((h_r - t_f_u) + y)/2. - y_cdg_r))/I_x_r

    I_x = np.pi*rad**3*t
    I_y = np.pi*rad**3*t
    I_xy = 0.

    theta = np.linspace(0,2*np.pi,101)


    data = []
    c = 0

    #Set origin of coordinates for the left-most figure (non-idealized)
    if w > max([h_l,h_r]):
        origin_x = -(max([h_l,h_r])/w*(max([h_l,h_r]) - w) + w)
    else:
        origin_x = -1.*max([h_l,h_r])
    origin_y = 0.

    x_coord_o = origin_x + (rad + t/2.)*np.cos(theta)
    y_coord_o = origin_y + rad + (rad + t/2.)*np.sin(theta)

    x_coord_i = origin_x + (rad - t/2.)*np.cos(theta)
    y_coord_i = origin_y + rad + (rad - t/2.)*np.sin(theta)  


    q_rad = lambda s: V_y/(np.pi*rad)*(np.cos(s) - 1)

    def x_sc_func(p):
        chi = p
        eq1 = np.sum(np.multiply(np.multiply(q_rad(theta),np.cos(theta)),[chi + x for x in rad*(np.cos(theta)+1)/2.]))
        return (eq1)
    # x_sc = integrate.quad(q_rad,0.,h_r)[0]*w/V_x

    x_sc = opt.fsolve(x_sc_func,[1.],args=())[0]

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1 = go.Scatter(
        x = x_coord_o,
        y = y_coord_o,
        mode='lines',
        name=r'Circ. Beam <br> V = {0:.2f} N <br> rad = {1:.0f} mm <br> t5 = {2:.0f} mm <br> A = {3:.0f} mm2 <br> I_x = {4:.0f} mm4 <br> I_y = {5:.0f} mm4 <br> I_xy = {6:.0f} mm4'.format(V_y,rad,t,A_tot,I_y,I_x,I_xy),
        hoverinfo='text',
        hoveron='fills',
        fill='toself', 
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    #Entitiy that plots the left-most (non-idealized) cross section
    trace1_i = go.Scatter(
        x = x_coord_i,
        y = y_coord_i,
        mode='lines',
        # hoveron='fills',
        hoverinfo='skip',
        fillcolor='rgba(26,150,65,0.)',
        line=dict(color='black'),
        visible=[True if c is 0 else False][0])

    #Entity that plots the location of the shear load application (shear center)
    trace1_sc = go.Scatter(
        x = [origin_x - rad + x_sc],
        y = [origin_y + rad],
        name = 'Shear Center',
        hovertext = r'Shear Center <br> x_sc = {0:.2f} <br> y_sc = {1:.2f} <br> w.r.t. exterior'.format(np.abs(x_sc),0.),
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'blue', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the location of the center of gravity
    trace1_cg = go.Scatter(
        x = [origin_x],
        y = [origin_y + rad],
        hovertext = r'Center of Gravity <br> x_cog = {0:.2f} <br> y_cog = {1:.2f} <br> w.r.t. exterior'.format(0.,0.),
        name = 'Center of Gravity',
        hoverinfo = 'text',
        mode = 'markers',
        marker = dict(color = 'purple', size = 10),
        visible=[True if c is 0 else False][0])
    
    #Entity that plots the neutral axis
    trace1_na = go.Scatter(
        x = np.linspace(origin_x - rad/2.,origin_x+ rad/2.,21),
        y = np.linspace(origin_y + y_cdg, origin_y + y_cdg, 21),
        name = 'Neutral Axis',
        mode='lines',
        hoveron = 'points',
        hoverinfo = 'name',
        line = dict(
            color = ('purple'),
            dash = 'dashdot'))

    #Determination of the shear flow fields  
    rad_par_flow_f = q_rad(theta)

    y_par_coord_f = [y_coord_o[x] - rad_par_flow_f[x]*np.sin(theta[x])*100*rad/V_y for x in range(len(rad_par_flow_f))]
    x_par_coord_f = [x_coord_o[x] - rad_par_flow_f[x]*np.cos(theta[x])*100*rad/V_y for x in range(len(rad_par_flow_f))]

    #Entity that plots the shear flow around the circle
    trace2_l = go.Scatter(
        x= x_par_coord_f,
        y= y_par_coord_f,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in rad_par_flow_f],
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    trace2_f = go.Scatter(
        x= x_par_coord_f + x_coord_o.tolist()[::-1],
        y= y_par_coord_f + y_coord_o.tolist()[::-1],
        mode='lines',
        fill = 'toself',
        hoverinfo='skip',
        # fillcolor = 'red',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    trace3 = go.Scatter(
        x = np.linspace(origin_x + rad - rad/4.,origin_x + rad + rad/4.,4),
        y = np.linspace(origin_y + rad, origin_y + rad, 4),
        name = 'Section Cut',
        hoverinfo = 'name',
        line = dict(
            color = ('black'),
            dash = 'dashdot'))

    #Set the origin of coordinates for the right-most figure (idealized)
    if w > max([h_l,h_r]):
        origin_xd = (max([h_l,h_r])/w*(max([h_l,h_r]) - w) + w)
    else:
        origin_xd = 1.*max([h_l,h_r])

    origin_yd = rad

    theta_d = np.linspace(0,2*np.pi,num+1).tolist()

    #Define the x and y coordinates of the right-most (idealized) cross section. The cross section is lumped to its midplane
    x_fu = origin_xd + rad*np.cos(theta_d)

    y_fu = origin_yd + rad*np.sin(theta_d)

    #Define x and y coordinates in a format so the discretized shear flow can be properly plotted
    x_fu = list(itertools.chain.from_iterable(itertools.repeat(x_fu[x], 2) for x in range(len(x_fu))))
    y_fu = list(itertools.chain.from_iterable(itertools.repeat(y_fu[x], 2) for x in range(len(y_fu))))

    x_coord_d = list(itertools.chain.from_iterable(itertools.repeat(x_coord_o[x], 2) for x in range(len(x_coord_o))))
    y_coord_d = list(itertools.chain.from_iterable(itertools.repeat(y_coord_o[x], 2) for x in range(len(y_coord_o))))

    rad_int = list(itertools.chain.from_iterable(itertools.repeat(theta_d[x], 2) for x in range(len(theta_d))))

    y_flange_u = q_rad

    #Define discretized shear flow by averaging calculated flow on non-idealized cross section
    q_fu_s =  [integrate.quad(y_flange_u,rad_int[2*x],rad_int[2*x+2])[0]/(rad_int[2*x + 2] - rad_int[2*x]) for x in range(int(len(rad_int)/2 -1))]
    q_fu = [0.] + [q_fu_s[int(x/2)] for x in range(int(len(x_fu)-2))] + [0.]

    q_fu_pl = [x*rad for x in q_fu]

    q_fu[0] = q_fu_s[0]
    q_fu[-1] = q_fu_s[-1]

    #Lines of code that smooth out the sharp peaks resulting from the direct discretization of the shear flow
    theta_d = list(itertools.chain.from_iterable(itertools.repeat(theta_d[x], 2) for x in range(len(theta_d))))

    side_angle = [np.arccos((x_fu[2*x+2] - x_fu[2*x+1])/(math.sqrt((x_fu[2*x+2] - x_fu[2*x+1])**2 + (y_fu[2*x+2] - y_fu[2*x+1])**2)))*180/np.pi for x in range(int(len(x_fu)/2 -1))]

    side_angle = [side_angle[x +1] for x in range(int(len(side_angle)/2))] + [side_angle[int(len(side_angle)/2) + x] for x in range(int(len(side_angle)/2))]

    side_angle = [0.] + side_angle + [0.]

    sign_angle = [np.sign(y_fu[2*x+2] - y_fu[2*x+1]) for x in range(int(len(x_fu)/2 -1))]
    sign_angle = [sign_angle[x +1] for x in range(int(len(sign_angle)/2))] + [-sign_angle[int(len(sign_angle)/2) + x] for x in range(int(len(sign_angle)/2))]
    sign_angle = [0.] + sign_angle + [0.]


    y_par_coord_d = [y_fu[x] - q_fu[x]*np.sin(theta_d[x])*100*rad/V_y  for x in range(len(q_fu) - 1)]
    x_par_coord_d = [x_fu[x] - q_fu[x]*np.cos(theta_d[x])*100*rad/V_y  for x in range(len(q_fu) - 1)]

    x_par_coord_dx = [np.abs(q_fu[2*x + 1] - q_fu[2*x])*100*rad/V_y*np.sin(2*np.pi/num)*(-np.abs(np.cos(side_angle[x]*np.pi/180.)))  for x in range(int(len(x_fu)/2))]
    y_par_coord_dx = [np.abs(q_fu[2*x + 1] - q_fu[2*x])*100*rad/V_y*np.sin(2*np.pi/num)*np.sin(side_angle[x]*np.pi/180.)*sign_angle[x]  for x in range(int(len(x_fu)/2))]


    for x in range(int(len(x_fu)/4)):
        x_par_coord_d[2*x+1] = x_par_coord_d[2*x+1] + x_par_coord_dx[x]
        y_par_coord_d[2*x+1] = y_par_coord_d[2*x+1] + y_par_coord_dx[x]
    for x in range(int(len(x_fu)/4) -1):
        x_par_coord_d[int(len(x_par_coord_d)/2) + 2*x + 1] = x_par_coord_d[int(len(x_par_coord_d)/2) + 2*x + 1] + x_par_coord_dx[int(len(x_par_coord_dx)/2) + x]
        y_par_coord_d[int(len(x_par_coord_d)/2) + 2*x + 1] = y_par_coord_d[int(len(x_par_coord_d)/2) + 2*x + 1] + y_par_coord_dx[int(len(x_par_coord_dx)/2) + x]

    #Calculate the corresponding boom areas based on the ratios of shear flow between two adjacent booms
    A_rad = np.zeros(num+1)
    
    theta_d_A = np.linspace(0,2*np.pi,num+1).tolist()

    for x in range(0,num):
        A_rad[x] += t*(2*np.pi*rad)/(num -1)/6.*(2 + q_rad(theta_d_A[x+1])/min([-1e-05,q_rad(theta_d_A[x])]))
        A_rad[x +1] += t*(2*np.pi*rad)/(num -1)/6.*(2 + q_rad(theta_d_A[x])/min([-1e-05,q_rad(theta_d_A[x+1])]))
        
    A_rad = A_rad.tolist()

    res_A = A_rad
    res_A = list(itertools.chain.from_iterable(itertools.repeat(res_A[x], 2) for x in range(len(res_A))))


    #Entitiy that plots the right-most (idealized) cross section
    trace1b = go.Scatter(
        x =  x_fu,
        y =  y_fu,
        mode='lines+markers',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['Ab = {0:.0f} mm2'.format(k) for k in res_A],
        line=dict(color='grey'),
        visible = [True if c is 0 else False][0])

    #Entity that plots the discretized shear flow on the upper flange
    trace2b = go.Scatter(
        x= x_par_coord_d,
        y= y_par_coord_d,
        mode='lines',
        name='Shear flow [N/mm]',
        hoverinfo = 'text',
        hoveron = 'points',
        hovertext = ['{0:.2f} N/mm'.format(np.abs(k)) for k in q_fu],
        #hovertemplate = '%{hovertext:.2f}',
        # fill = 'toself',
        line=dict(color='red'),
        visible = [True if c is 0 else False][0])

    trace2b_f = go.Scatter(
        x= x_fu + x_par_coord_d[::-1],
        y= y_fu + y_par_coord_d[::-1],
        mode='lines',
        fill = 'toself',
        hoveron = 'fills',
        hoverinfo='skip',
        # fillcolor = 'red',
        line=dict(color='red'),
        visible=[True if c is 0 else False][0])

    trace3b = go.Scatter(
        x = np.linspace(origin_xd + rad - rad/4.,origin_xd + rad + rad/4.,4),
        y = np.linspace(origin_yd, origin_yd, 4),
        name = 'Section Cut',
        hoverinfo = 'name',
        line = dict(
            color = ('black'),
            dash = 'dashdot'))


    #Add all the plotting entities to a data list that maps the plotting
    data_temp = [trace1_i, trace2_f, trace2_l, trace1,trace1_sc, trace1_cg,trace1_na, trace3,trace2b, trace2b_f, trace1b, trace3b]

    data = data + data_temp
    c = c +1

    theta_arr = np.linspace(0,2*np.pi,20)

    annotation_list = []
    for x in theta_arr:
        annotation_list_temp = dict(x = origin_x + rad*np.cos(x) ,y = origin_y + rad + rad*np.sin(x),xref = "x", yref = "y", showarrow = True, arrowhead=2, arrowsize=1, arrowwidth=1, arrowcolor='red', ax= 15.*np.sin(x), ay = 15.*np.cos(x))
        annotation_list.append(annotation_list_temp)


    #Define layout options of figure
    layout = go.Layout(
        showlegend=False,
        hovermode= 'closest',
        #Add annotations with arrows indicating the shear flow direction and the shear load application which are proportional to the cross-section size
        annotations = annotation_list + 
                     
                       [dict(x = origin_x - rad + x_sc, y = origin_y + rad + 0.05*rad*2, xref = "x", yref = "y", text = "V",  font = dict( color = "black",  size = 14 ),showarrow = True, arrowhead=1, arrowsize=1.5, arrowwidth=1, arrowcolor='black', ax=0., ay= -40, xanchor = "center" ) ],
      #Define axes visibility properties
        xaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            showticklabels=False
        ),
        yaxis=dict(
            autorange=True,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks='',
            scaleanchor="x", 
            scaleratio=1,
            showticklabels=False
        )
    )

    return (data, layout)

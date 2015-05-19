'''
Integration of a reaction-diffusion system exhibiting anomalous diffusion with a power-law memory kernel

This file initializes a set of parameters, runs the integration for a relatively generous amount of time, and then plots the output

William Gilpin, Spakowitz Group at Stanford University, 2015
'''

## pick parameter values
params = dict()
# width of reactive well
params['WELL_DIAM'] = .2
# width of overall potential well. smaller this is, the stronger the forcing
params['POT_DIAM'] = .05
params['KAPPA'] = 4e-1
params['DCOEFF'] = 1e-4
# alpha less than one
params['ALPHA']= .5


## set integrator settings
settings = dict()

space_pts = 25
ACTUAL_LENGTH = .5
dx = ACTUAL_LENGTH/space_pts
space = linspace(0.0, ACTUAL_LENGTH, space_pts)
space = space+dx

time_pts = 1e5
start_time = 0.0
stop_time = 17.0
dt = (stop_time-start_time)/time_pts
times = linspace(start_time, stop_time, time_pts)
times = times + dt

# initial conditions
y0 = ones(space_pts)

# settings['dx'] = dx
settings['space'] = space
settings['times'] = times


sol = rk4_mesh(y0, nxt_step, settings, params)


# Plot the results

# plot space-time concentration
figure()
imshow(sol, aspect='auto')

# plot time slices
figure()
solrange = len(sol[0,:])
num_slices = 50
hold(True)
slice_range = floor(solrange/num_slices)
plot(space, sol[:,::slice_range])

plot(space, sol[:,-1])
xlabel("Position along well")
ylabel("Probability")
# 'alpha_' + str(params['ALPHA']) + '__diff_' + str(params['DCOEFF'])
# savefig('harmonic_alpha_1__a_p2__V_p4.pdf')


# figure()
total_conc = sum(sol, axis = 0)
plot(times,total_conc[:-1])
ylabel("Total count")
xlabel("time")
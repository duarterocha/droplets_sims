from scripts import *

class DropletTempProblemMafromRa(DropletTempProblem): # Inherit from the previous problem

    def __init__(self):
        super(DropletTempProblemMafromRa, self).__init__()

        # default 2% of the cartesian 2d are covered by positive stream function, 98% negative values
        self.positive_stream_volume_fraction=0.02

        self._initial_Ma=0 # Initial Ma value, must be stored when set_Ma is called before the problem is initialized

    # include getter and setter for Ma unknown

    def get_Ma(self, symbolic=False):
        if symbolic:
            return var("Ma", domain="globals")  # symbolic is just the variable value of the Lagrange multiplier Ma
        elif not self.is_initialised():
            return self._initial_Ma  # when the problem is not initialized, there is no degree of freedom Ra yet
        else:
            return self.get_ode("globals").get_value("Ma", as_float=True)  # Get the current value

    def set_Ma(self,value):
        if self.is_initialised():
            self.get_ode("globals").set_value(Ma=value) # set the current value
        else:
            self._initial_Ma=value # set the initial value

    def actions_after_newton_step(self):
        super(DropletTempProblemMafromRa, self).actions_after_newton_step()
        pos_area=self.get_mesh("droplet").evaluate_observable("positive_stream_area_fraction")
        print("Currently at Ma=",self.get_Ma(),"Ra=",self.Ra.value,"POSITIVE COVERED AREA",pos_area,"DESIRED",self.positive_stream_volume_fraction)
        # Set the Ma parameter (which is unused for the solution here) to the current value
        # it is necessary for the output
        self.Ma.value=self.get_Ma()

    def define_problem(self):
        super(DropletTempProblemMafromRa, self).define_problem()
        # Ma is now not a parameter anymore, but an unknown
        # First of all, we must define this unknown. To that end, we define a global Lagrange multiplier
        # There is nothing (i.e. 0) added to the residual
        MaDef=GlobalLagrangeMultiplier(Ma=0)
        MaDef+=InitialCondition(Ma=self._initial_Ma) # Set the initially selected value
        MaDef+=DirichletBC(Ma=self._initial_Ma) # Fix the initial value (used in the first solve)
        self.add_equations(MaDef@"globals") # and add it to a domain "globals", which is an ODE domain
        # We must define an equation for Ma to determine it. As described above, it is a Lagrange multiplier to ensure
        #   integral_2d_cartesian(heaviside(streamfunc)-phi)=0
        # However, since
        #   heaviside(x)=(1 if x>0 else 0)
        # is not differentiable, it is better to smooth the slope between 0 and 1 a bit:
        # TODO: It would be beneficial to adjust he large factor depending on the range of the stream function
        # If you make the factor too small, it won't work well for small Ma,Ra
        # Too large factors may hamper the solution at higher Ma,Ra
        smooth_heaviside=lambda x : (tanh(x*10000000.0)+1)/2
        # add the condition in a Lagrange multiplier sense,
        # so that the appropriate Rayleigh number is indeed found if
        #   integral_2d_cartesian(heaviside(streamfunc)-phi)=0
        MaEq=WeakContribution(smooth_heaviside(var("streamfunc"))-self.positive_stream_volume_fraction,testfunction("Ma",domain="globals"),coordinate_system=cartesian)@"droplet"
        self.add_equations(MaEq)

with DropletTempProblemMafromRa() as problem:
    problem.set_c_compiler("tcc")

    # Select the contact angle
    problem.contact_angle=60*degree
    problem.resolution=0.025 # For the final publication, make it fine (i.e. small)

    # Start with a small Ra
    problem.Ra.value=1

    # And set a reasonable starting Ma
    # It is always good to start in the Ma-vs-Ra regime. There, the area with positive stream function changes with Ma
    # so that the solver can converge. If you are in pure Ra or pure Ma domain, the positive stream function are won't
    # change with Ra and the solver has a hard time
    problem.positive_stream_volume_fraction = 0.01  # One percent of the area covered by positive stream function
    problem.set_Ma(20)

    # First solve. Ra is pinned by the DirichletBC to the initial value of Ra.
    # Thereby we solve for the flow, etc, but do not adjust Ra yet
    problem.solve()

    # Now when all other fields are converged, let's activate the Ra finding
    problem.get_ode("globals").set_dirichlet_active(Ma=False)  # Remove the fixation of Ma, now Ma will be adjusted
    problem.reapply_boundary_conditions()  # Since the boundary conditions are changed (Ma is now free), we must call this

    # The first solve is hard. Let it converge globally by line search
    # and allow for many iterations. It will only work if a reasonable initial guess in terms of Ma,Ra is provided
    problem.solve(globally_convergent_newton=True, max_newton_iterations=100)
    # We found Ma(Ra) that fulfills the condition, output!
    problem.output()

    # Create a file in the output directory to write the curve Ra(Ma)
    parameter_curve_file = open(problem.get_output_directory(relative_path="parameter_curve.txt"), "w")
    parameter_curve_file.write("#Ma\tRa\n")  # write header
    # write Ma, Ra and flush to write instantly (i.e. do not buffer)
    parameter_curve_file.write("\t".join(map(str, [problem.get_Ma(), problem.Ra.value])) + "\n")
    parameter_curve_file.flush()

    # Remove the gas phase from solving (it is always the same!)
    with problem.select_dofs() as dofs:
        dofs.unselect("gas")  # Gas is not solved for

        # Perform a gradual increase of Ra, always with the constraint that Ra is adjusted according to the condition
        ds = 0.1  # Initial step (first step is the step in Ra, all other steps are in terms of the arc length)
        while problem.Ra.value < 100000:  # Limit Ma
            ds = problem.arclength_continuation("Ra", ds)  # Perform an adaptive arclength continuation

            # Optionally. If set, we can leave the current branch of the solution
            problem.reset_arc_length_parameters()  # and ds will be again in terms of Ra, not the arclength

            problem.output_at_increased_time()  # Output and write to the parameter curve file
            parameter_curve_file.write("\t".join(map(str, [problem.get_Ma(), problem.Ra.value])) + "\n")
            parameter_curve_file.flush()
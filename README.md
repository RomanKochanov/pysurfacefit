<h1>Multidimensional data fitting tool in Python</h1>

Python package for linear and non-linear fitting of multidimensional data.<br>

PySurfaceFit was initially planned as a narrow-scope tool for fitting the molecular potential energy surface (PES) models to the <it>ab initio</it> data computed via quantum chemistry methods. Currently, PySurfaceFit is a flexible general-purpose least squares fitting tool supporting a number of global and local optimization methods.<br>

PySurfaceFit supports fit constraints and <a href="https://en.wikipedia.org/wiki/Ridge_regression#Tikhonov_regularization">Tikhonov regularization</a> for the data points and models parameters. For the Sympy models, based on the <a href="https://www.sympy.org">Sympy</a> symbolic math library, is is possible to use analytic Jacobians and Hessians if required by the particular optimization method.

Currently, PySurfaceFit heavily relies on the  <a href="https://docs.scipy.org/doc/scipy/reference/optimize.html">Scipy optimization and root finding library</a> (scipy.optimize).<br>

<b>The basic supported <it>local</it> optimization methods:</b>

<ol>
  <li>Levenberg-Marquardt</li>
  <li>Trust Region Reflective</li>
  <li>Nelder-Mead</li>
  <li>Powell</li>
  <li>CG (Conjugate Gradient)</li>
  <li>Newton-CG</li>
  <li>Broyden-Fletcher-Goldfarb-Shanno (BFGS)</li>
  <li>L-BFGS-B</li>
</ol>

<b>The basic supported <it>global</it> optimization methods:</b>

<ol>
  <li>Simulated Annealing</li>
  <li>Basin Hopping Optimization</li>
</ol>

Other custom optimization method backends are fairly easy to implement (see the Notebook tutorials on the optimization methods integration).

<h1>Installation</h1>

<h2>Pulling from Github</h2>

<ol>
  <li>git clone https://github.com/RomanKochanov/pysurfacefit.git</li>
  <li>cd pysurfacefit</li>
  <li>pip install .</li>
</ol> 

<h2>From Python Package index</h2>

pip install pysurfacefit

<h1>PySurfaceFit CLI tool</h1>

PySurfaceFit provides two options of usage: as a CLI tool, and as a Python library. The CLI tool provides almost all basic functionality for managing the fit data points, creating fit models, and displaying the results of a fit. For beginners, it is recommended using the CLI tool to avoid most of the programming fuss. 

<b>Once installed, PySurfaceFit can be run from the console:</b>

<pre><code>
pysurfacefit
</code></pre>

<b>To obtain initial help, use the --help flag:</b>

<pre><code>
pysurfacefit --help
</code></pre>

<b>OUTPUT:</b>

<pre>usage: pysurfacefit [-h] [--config CONFIG] [--startproject STARTPROJECT]
                    [--template TEMPLATE [TEMPLATE ...]] [--merge MERGE]
                    [--set SET] [--init] [--split] [--fit] [--stat]
                    [--codegen] [--plot [PLOT [PLOT ...]]] [--calc]

Python package for linear and non-linear fitting of multidimensional data.

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG       Configuration file (mandatory for all steps except
                        --startproject)
  --startproject STARTPROJECT
                        Stage 1: start empty project, create dummy config
  --template TEMPLATE [TEMPLATE ...]
                        _________1a: us a template for new project
  --merge MERGE         _________1b: take defaults from a config file
  --set SET             _________1c: explicitly set project parameters
  --init                Stage 2: create model package
  --split               Stage 3: split initial datafile with weighting scheme
  --fit                 Stage 4: start fitting model to data
  --stat                Stage 5: calculate statistics
  --codegen             Stage 6: generate Fortran code for fitted model
  --plot [PLOT [PLOT ...]]
                        Stage 7: plot sections of the model and compare to
                        data
  --calc                Stage 8: calculate model values on grid
</pre>

<b>This will generate an emply project with the sample config file, which will look something like this:</b><br>

<pre><font color="#34E2E2">####################</font>
<font color="#34E2E2"># GENERAL SETTINGS #</font>
<font color="#34E2E2">####################</font>
<font color="#FFD7D7">[GENERAL]</font>

<font color="#34E2E2"># Name of the project.</font>
project: myproj

<font color="#34E2E2">##############################</font>
<font color="#34E2E2"># FITTING DATA SPECIFICATION #</font>
<font color="#34E2E2">##############################</font>
<font color="#FFD7D7">[DATA]</font>

<font color="#34E2E2"># CSV file containing fitting data.</font>
datafile:

<font color="#34E2E2"># Define rules to split the datafile to dataspec.</font>
<font color="#34E2E2"># Should be ignored if dataspec is empty.</font>
split_column:
split_values:
split_weights:

<font color="#34E2E2"># Resulting data specification split-file.</font>
<font color="#34E2E2"># N.B.: if empty, then the raw datafile is used.</font>
<font color="#34E2E2"># Data specification *.txt file has the following format:</font>
<font color="#34E2E2"># alias   path       wht_mul   type    include</font>
<font color="#34E2E2"># ------- ---------  --------  ------  -------</font>
<font color="#34E2E2"># test1   test1.csv  1         0       1       </font>
<font color="#34E2E2"># test2   test2.csv  1         0       1</font>
<font color="#34E2E2"># test3   test3.csv  1         0       1</font>
dataspec: dataspec.txt

<font color="#34E2E2"># Names and units of the input data columns. </font>
<font color="#34E2E2"># Input names must be separated by semicolon </font>
<font color="#34E2E2"># and should not contain dots.</font>
<font color="#34E2E2"># E.g.: X;Y;Z for names, X:X_UNIT;Y:Y_UNIT for units</font>
input_columns:
input_units:

<font color="#34E2E2"># Name of the output data column. See explanation above.</font>
output_column:
output_units:

<font color="#34E2E2"># Weight function.</font>
wht_fun: lambda v: <font color="#AD7FA8">1</font>

<font color="#34E2E2"># Global data cutoff</font>
global_cutoff_max:
global_cutoff_min:

<font color="#34E2E2">################################</font>
<font color="#34E2E2"># FITTING MODEL SPECIFICATIONS #</font>
<font color="#34E2E2">################################</font>
<font color="#FFD7D7">[MODEL]</font>

<font color="#34E2E2"># Model package name.</font>
model: fitmodel

<font color="#34E2E2"># Model arguments. Argument names must be in the same order</font>
<font color="#34E2E2"># as the data column names from the DATA section.</font>
<font color="#34E2E2"># E.g.: X;Y;Z</font>
arguments:

<font color="#34E2E2">##################</font>
<font color="#34E2E2"># FITTER OPTIONS #</font>
<font color="#34E2E2">##################</font>
<font color="#FFD7D7">[FIT]</font>

<font color="#34E2E2"># Weighted / unweighted fit mode.</font>
weighted_fit: True

<font color="#34E2E2"># Use &quot;rubber&quot; regularization.</font>
rubber_on: True

<font color="#34E2E2"># Fitting method (trf, lm, basinhopping).</font>
<font color="#34E2E2"># Available &quot;local&quot; methods:</font>
<font color="#34E2E2">#    &quot;trf&quot; -&gt; Trust Region Reflective (least squares, bounded).</font>
<font color="#34E2E2">#    &quot;lm&quot;  -&gt; Levenberg Marquardt (least squares, unbounded).</font>
<font color="#34E2E2">#    &apos;Nelder-Mead&apos;  -&gt; Gradient-free descent Nelder-Mead method.</font>
<font color="#34E2E2">#    &apos;Powell&apos;  -&gt; modified Powell algorithm.</font>
<font color="#34E2E2">#    &apos;CG&apos;  -&gt; conjugate gradient algorithm.</font>
<font color="#34E2E2">#    &apos;BFGS&apos;  -&gt; BFGS algorithm.</font>
<font color="#34E2E2">#    &apos;Newton-CG&apos;  -&gt; Newton-CG algorithm.</font>
<font color="#34E2E2">#    &apos;L-BFGS-B&apos;  -&gt; L-BFGS-B algorithm.</font>
<font color="#34E2E2">#    &apos;TNC&apos;  -&gt; truncated Newton (TNC) algorithm.</font>
<font color="#34E2E2">#    &apos;COBYLA&apos;  -&gt; Constrained Optimization BY Linear Approximation.</font>
<font color="#34E2E2">#    &apos;SLSQP&apos;  -&gt; Sequential Least Squares Programming.</font>
<font color="#34E2E2">#    &apos;trust-constr&apos;  -&gt; minimize a scalar function subject to constraints.</font>
<font color="#34E2E2">#    &apos;dogleg&apos;  -&gt; dog-leg trust-region algorithm.</font>
<font color="#34E2E2">#    &apos;trust-ncg&apos;  -&gt; Newton conjugate gradient trust-region algorithm.</font>
<font color="#34E2E2">#    &apos;trust-exact&apos;  -&gt; &quot;nearly exact trust-region algorithm&quot;</font>
<font color="#34E2E2">#    &apos;trust-krylov&apos;  -&gt; &quot;early exact trust-region algorithm&quot;</font>
<font color="#34E2E2"># Available &quot;global&quot; methods:</font>
<font color="#34E2E2">#    &quot;basinhopping&quot;  -&gt; local descents from a random starting point (bounded)</font>
<font color="#34E2E2">#    &quot;anneal&quot; -&gt; Simulated Annealing method</font>
fitting_method: trf

<font color="#34E2E2"># Calculate analytic jacobian (only if ModelSympy is used).</font>
analytic_jacobian: True

<font color="#34E2E2"># Fit options. Must be given as a list separated by semicolon.</font>
<font color="#34E2E2"># The valid options for most of the supported fitting methods</font>
<font color="#34E2E2"># can be found in the corresponding scipy.optimize documentation.  </font>
<font color="#87FFAF">fit_options: max_nfev=</font><font color="#AD7FA8">200</font>;

<font color="#34E2E2">##################################</font>
<font color="#34E2E2"># STATISTICS CALCULATION OPTIONS #</font>
<font color="#34E2E2">##################################</font>
<font color="#FFD7D7">[STAT]</font>

<font color="#34E2E2"># Fit statistics file.</font>
stat_file: stat.out

<font color="#34E2E2"># Calculate outlier statistics.</font>
outlier_stats_flag: True

<font color="#34E2E2"># Type of statistics: cook, dffits, leverage ,student.</font>
outlier_stats_type: cook

<font color="#34E2E2"># Output generated symbolic function, for debug purposes only.</font>
output_symbolic_func:

<font color="#34E2E2">##########################</font>
<font color="#34E2E2"># CODE GENERATOR OPTIONS #</font>
<font color="#34E2E2">##########################</font>
<font color="#FFD7D7">[CODEGEN]</font>

<font color="#34E2E2"># Create Fortran code.</font>
create_fortran: False

<font color="#34E2E2"># Compare original Python code with the generated Fortran.</font>
compare_fortran: False

<font color="#34E2E2"># Fortran compiler executable</font>
compiler_fortran: ifort

<font color="#34E2E2"># Grid specifications to calculate on.</font>
<font color="#34E2E2"># Format of the specification must be as follows: </font>
<font color="#34E2E2"># X=XMIN:XSTEP:XMAX; Y=YMIN:YSTEP:YMAX; ...</font>
gridspec:

<font color="#34E2E2">####################</font>
<font color="#34E2E2"># PLOTTING OPTIONS #</font>
<font color="#34E2E2">####################</font>
<font color="#FFD7D7">[PLOTTING]</font>

<font color="#34E2E2"># ATTENTION: names of the plot coordinates correspond to the </font>
<font color="#34E2E2"># model argumen names given in the MODEL section.</font>

<font color="#34E2E2"># Type of the plot.</font>
<font color="#34E2E2"># Available plot modes:</font>
<font color="#34E2E2">#   =&gt; &quot;residuals&quot;: </font>
<font color="#34E2E2">#           plot unweighted fit residuals.</font>
<font color="#34E2E2">#   =&gt; &quot;sections&quot;: </font>
<font color="#34E2E2">#           plot sections/cuts of the fitted model vs datapoints.</font>
plot_mode: residuals

<font color="#34E2E2"># Plot coordinate grid specifications.</font>
<font color="#34E2E2"># Format of the specification must be as follows: </font>
<font color="#34E2E2"># X=XVAL; Y=YMIN:YSTEP:YMAX; ...</font>
<font color="#34E2E2"># In case of the fixed coordinate, variable would only have a value, e.g. X=XVAL.</font>
<font color="#34E2E2"># In case of unfixed coordinate, there should be a full grid, e.g. Y=YMIN:YSTEP:YMAX</font>
<font color="#34E2E2"># N.B.: 1) order of the unfixed coords affects the order of the plot axes.</font>
<font color="#34E2E2">#       2) order of the binding MUST correspond to the argument of the model&apos;s </font>
<font color="#34E2E2">#          __func__ method.</font>
gridspec:

<font color="#34E2E2"># Model components to plot.</font>
model_components:

<font color="#34E2E2"># Calculate model&apos;s components.</font>
calculate_components: False

<font color="#34E2E2"># Plot outlier statistics in color. If False, each datagroup has its own color.</font>
plot_outlier_stats: False

<font color="#34E2E2"># Plot weighted residuals.</font>
resids_weighted: False

<font color="#34E2E2"># X axes to plot the resuduals versus. </font>
<font color="#34E2E2"># If empty, defaults to [DATA][output_column].</font>
resids_x_axes:

<font color="#34E2E2"># Scatter settings (2D, 3D case)</font>
scatter_opacity: <font color="#AD7FA8">1</font>
marker_size: <font color="#AD7FA8">20</font>
resize_by_weights: False

<font color="#34E2E2"># Surface settings (3D case)</font>
surface_opacity: <font color="#AD7FA8">0.4</font>

<font color="#34E2E2">#######################</font>
<font color="#34E2E2"># CALCULATION OPTIONS #</font>
<font color="#34E2E2">#######################</font>
<font color="#FFD7D7">[CALC]</font>

<font color="#34E2E2"># Output file to output calculated model. </font>
<font color="#34E2E2"># Column names are defined in DATA section.</font>
output_file: calc.csv

<font color="#34E2E2"># Grid specifications to calculate on.</font>
<font color="#34E2E2"># Format of the specification must be as follows: </font>
<font color="#34E2E2"># X=XMIN:XSTEP:XMAX; Y=YMIN:YSTEP:YMAX; ...</font>
gridspec:
</pre>

<b>Now, let's perform a test fitting case step by step. First, let's generate a three-dimensional data we will fit our model to.</b>


```python
import numpy as np
from pysurfacefit.grids import Grid

# set 1D grids for each coordinate
x = np.arange(-10,10,1.0)
y = np.arange(-15,-5,0.1)
z = np.arange(0,20,5.0)

# Create a 3D grid
grid = Grid('x',x,'y',y,'z',z)

# Define a test function
func = lambda x,y,z: x**2 + 0.1*y**3 + 10*z**4 + 10

# Calculate a function on the 3D grid
v = grid.calculate(func).flatten()

# Save the resulting data to CSV file
x, y, z = grid.get_meshes(flat=True)
np.savetxt('sampledata.csv',list(zip(x,y,z,v)),delimiter=';',header='x;y;z;v',comments='')
```

<br>

<b>After the test data file is generated, we can create an empty project named "test" from it.</b><br>

<pre>pysurfacefit --startproject test</pre><br>

<b>To generate a ready-to-go configuration file, one must provide additional data such as names for the argument columns (or "inputs"), and the name of the function column (or "output"). The entered names must correspond to the column names in the data file:</b><br>

<pre>Enter the name of the data points CSV file: sampledata.csv
Enter semicolon-separated names for model inputs: x; y; z 
Enter name for model output: v
</pre><br>

<b>If everything went fine, we should have the new fitting project created in a nested subfolder.</b><br>

<pre>Created new project: test
New config file has been added: test/config.ini
Sample data specification file has been added: dataspec.txt
</pre><br>

<b>Now the project folder "test" inlcudes only three files:</b><br>

<pre>cd test; ls
config.ini  dataspec.txt  sampledata.csv</pre><br>

<b>The file "config.ini" is the main configuration file containing the project settings (see the example of the file above). This file is separated into the following sections:</b>

<ol>
  <li><b>GENERAL</b></li>
    <em>Basic settings: contains only the name of the project so far.</em><br><br>
  <li><b>DATA</b></li>
    <em>Settings for the data to fit: defines the source file for the fitted data, names and units of the input and output columns, weight function and cutoff. Also contains settings for data splitting to apply more advanced weighting schemes which are defined in the dataspec.txt file.</em><br><br>
  <li><b>MODEL</b></li>
    <em>Containes the model name and names for model arguments. Usually, the argument names are the same as the names of the data inputs defined in the DATA section.</em><br><br>
  <li><b>FIT</b></li>
    <em>Settings for the fitter. Containes the switches for turning on/off weights and regularization, name of the optimization method, switch for analytic Jacobian, and number of iterations setting.</em><br><br>
  <li><b>STAT</b></li>
    <em>Additional settings to display the fit statistics.</em><br><br>
  <li><b>CODEGEN</b></li>
    <em>Settings for the code generation from the fitted model. Allows choosing the compiler. Currently only the Fortan code generation is supported.</em><br><br>
  <li><b>PLOTTING</b></li>
    <em>Settings for plotting the fit residuals and 1D/2D sections of the fitted model. 3D volume visualization is planned to be included.</em><br><br>
  <li><b>CALC</b></li>
    <em>Defines the output file and grid to calculate the fited model on.</em><br><br>
</ol>

<b>To fit the data in the sampledata.csv file, let's generate a model file.</b>

<pre>
pysurfacefit --config config.ini --init
</pre>

Note, that for all tasks except --startproject, we must supply the configuraion file using the --config key.<br><br>

<b>The model is created and stored in a separate Python module file. The name of this file is defined in the MODEL section using the "model" parameter.</b>

<pre>
ls
config.ini  dataspec.txt  fitmodel.py  sampledata.csv
</pre>


<b>The default model has just one fitting parameter which is a constant. Here is what a default model file looks like:</b>

<pre><font color="#5FD7FF">import</font> os

<font color="#5FD7FF">from</font> pysurfacefit.models.sympy <font color="#5FD7FF">import</font> ModelSympy
<font color="#5FD7FF">from</font> pysurfacefit.fitpars <font color="#5FD7FF">import</font> Par, Parameters

<font color="#FCE94F">class</font> <font color="#34E2E2"><b>fitmodel</b></font>(ModelSympy):
    <font color="#AD7FA8">&quot;&quot;&quot;</font>
<font color="#AD7FA8">    Test model for fitting data with just one parameter.</font>
<font color="#AD7FA8">    &quot;&quot;&quot;</font>
    <font color="#FCE94F">def</font> <font color="#34E2E2"><b>__init__</b></font>(self,calc_switch=<font color="#AD7FA8">&apos;numbified&apos;</font>):
        self.__check_symbolic__ = <font color="#34E2E2"><b>False</b></font>
        self.__calc_switch__ = calc_switch <font color="#34E2E2"># symbolic, lambdified, numbified</font>
        self.__components__ = {}

        <font color="#34E2E2"># Initialize empty parameters.</font>
        self.__params__ = Parameters()

        <font color="#34E2E2"># Add constant parameter.</font>
        self.__params__.append(group=<font color="#AD7FA8">&apos;constant&apos;</font>, pars=[
            Par(name=<font color="#AD7FA8">&apos;constant&apos;</font>, value=<font color="#AD7FA8">0.0</font>, flag=<font color="#34E2E2"><b>True</b></font>),
        ])

    <font color="#FCE94F">def</font> <font color="#34E2E2"><b>__units__</b></font>(self):
        <font color="#FCE94F">return</font> {<font color="#AD7FA8">&quot;input&quot;</font>: {}, <font color="#AD7FA8">&quot;output&quot;</font>: <font color="#AD7FA8">&quot;None&quot;</font>}
        <font color="#FCE94F">def</font> <font color="#34E2E2"><b>__func__</b></font>(self,params,x,y,z):

        <font color="#34E2E2"># constant</font>
        constant = params[<font color="#AD7FA8">&apos;constant&apos;</font>].get_value()

        <font color="#34E2E2"># final result</font>
        res = constant

        <font color="#34E2E2"># components</font>
        self.__components__[<font color="#AD7FA8">&apos;constant&apos;</font>] = constant

        <font color="#FCE94F">return</font> res

model = fitmodel()

<font color="#FCE94F">if</font> os.path.exists(<font color="#AD7FA8">&apos;fitmodel.csv&apos;</font>):
    model.load_params(<font color="#AD7FA8">&apos;fitmodel.csv&apos;</font>)
<font color="#FCE94F">else</font>:
    model.save_params(<font color="#AD7FA8">&apos;fitmodel.csv&apos;</font>)
</pre>

<b>So far, we will use this default constant model to fit the three-dimensional data we have generated earlier. To do that, use the "fit" command:</b><br>

<pre>pysurfacefit --config config.ini --fit</pre><br>

<b>This will start the fitting process and will print the intermediate fit statistics to the STDOUT:</b><br>

<pre>
<font size="-2">
Treating sampledata fitgroup as FitPoints
BEGIN FIT
USING SCIPY.OPTIMIZE.LEAST_SQUARES: METHOD= trf
METHOD OPTIONS: {&apos;max_nfev&apos;: 200, &apos;bounds&apos;: [(-inf,), (inf,)], &apos;jac&apos;: &lt;function Fitter.fit_least_squares.&lt;locals&gt;.&lt;lambda&gt; at 0x7f84462e54d0&gt;}

==========================================
&lt;&lt;&lt;&lt; calling __sympy_initialize_func__ &gt;&gt;&gt;&gt;
==========================================
Progress:
     - creating sympy objects for inputs
     - creating sympy objects for parameters
     - get the Sympy expression by calling function with the Sympy objects
     - create lambdified Python function from sympy expression
     - create compiled (numbified) code from the lambdified function
==========================================

CALC FUN     1&gt;&gt;&gt;  DIST:0.000000000e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.324527156e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.325e+14   UNWHT_SSE: 5.325e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================

==========================================
&lt;&lt;&lt;&lt; calling __sympy_initialize_jac__ &gt;&gt;&gt;&gt;
==========================================
Progress:
     - get the Sympy expression by calling function with the Sympy objects
     - create lambdified Python function from sympy expression
     - create compiled (numbified) code from the lambdified function
==========================================

CALC JAC     1&gt;&gt;&gt;  DIST:0.000000000e+00
CALC FUN     2&gt;&gt;&gt;  DIST:1.000000000e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.324502669e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.325e+14   UNWHT_SSE: 5.325e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================
CALC JAC     2&gt;&gt;&gt;  DIST:1.000000000e+00
CALC FUN     3&gt;&gt;&gt;  DIST:3.000000000e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.324453696e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.324e+14   UNWHT_SSE: 5.324e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================
CALC JAC     3&gt;&gt;&gt;  DIST:3.000000000e+00
CALC FUN     4&gt;&gt;&gt;  DIST:7.000000000e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.324355753e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.324e+14   UNWHT_SSE: 5.324e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================
CALC JAC     4&gt;&gt;&gt;  DIST:7.000000000e+00
CALC FUN     5&gt;&gt;&gt;  DIST:1.500000000e+01  SSE_RUB:0.000000000e+00 SSE_TOT:5.324159873e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.324e+14   UNWHT_SSE: 5.324e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================
CALC JAC     5&gt;&gt;&gt;  DIST:1.500000000e+01
CALC FUN     6&gt;&gt;&gt;  DIST:3.100000000e+01  SSE_RUB:0.000000000e+00 SSE_TOT:5.323768145e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.324e+14   UNWHT_SSE: 5.324e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================
CALC JAC     6&gt;&gt;&gt;  DIST:3.100000000e+01
CALC FUN     7&gt;&gt;&gt;  DIST:6.300000000e+01  SSE_RUB:0.000000000e+00 SSE_TOT:5.322984811e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.323e+14   UNWHT_SSE: 5.323e+14   WHT_SD: 2.579e+05   UNWHT_SD: 2.579e+05
===============================
CALC JAC     7&gt;&gt;&gt;  DIST:6.300000000e+01
CALC FUN     8&gt;&gt;&gt;  DIST:1.270000000e+02  SSE_RUB:0.000000000e+00 SSE_TOT:5.321418635e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.321e+14   UNWHT_SSE: 5.321e+14   WHT_SD: 2.579e+05   UNWHT_SD: 2.579e+05
===============================
CALC JAC     8&gt;&gt;&gt;  DIST:1.270000000e+02
CALC FUN     9&gt;&gt;&gt;  DIST:2.550000000e+02  SSE_RUB:0.000000000e+00 SSE_TOT:5.318288249e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.318e+14   UNWHT_SSE: 5.318e+14   WHT_SD: 2.578e+05   UNWHT_SD: 2.578e+05
===============================
CALC JAC     9&gt;&gt;&gt;  DIST:2.550000000e+02
CALC FUN    10&gt;&gt;&gt;  DIST:5.110000000e+02  SSE_RUB:0.000000000e+00 SSE_TOT:5.312035342e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.312e+14   UNWHT_SSE: 5.312e+14   WHT_SD: 2.577e+05   UNWHT_SD: 2.577e+05
===============================
CALC JAC    10&gt;&gt;&gt;  DIST:5.110000000e+02
CALC FUN    11&gt;&gt;&gt;  DIST:1.023000000e+03  SSE_RUB:0.000000000e+00 SSE_TOT:5.299560985e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.300e+14   UNWHT_SSE: 5.300e+14   WHT_SD: 2.574e+05   UNWHT_SD: 2.574e+05
===============================
CALC JAC    11&gt;&gt;&gt;  DIST:1.023000000e+03
CALC FUN    12&gt;&gt;&gt;  DIST:2.047000000e+03  SSE_RUB:0.000000000e+00 SSE_TOT:5.274738099e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.275e+14   UNWHT_SSE: 5.275e+14   WHT_SD: 2.568e+05   UNWHT_SD: 2.568e+05
===============================
CALC JAC    12&gt;&gt;&gt;  DIST:2.047000000e+03
CALC FUN    13&gt;&gt;&gt;  DIST:4.095000000e+03  SSE_RUB:0.000000000e+00 SSE_TOT:5.225595644e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.226e+14   UNWHT_SSE: 5.226e+14   WHT_SD: 2.556e+05   UNWHT_SD: 2.556e+05
===============================
CALC JAC    13&gt;&gt;&gt;  DIST:4.095000000e+03
CALC FUN    14&gt;&gt;&gt;  DIST:8.191000000e+03  SSE_RUB:0.000000000e+00 SSE_TOT:5.129324001e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.129e+14   UNWHT_SSE: 5.129e+14   WHT_SD: 2.532e+05   UNWHT_SD: 2.532e+05
===============================
CALC JAC    14&gt;&gt;&gt;  DIST:8.191000000e+03
CALC FUN    15&gt;&gt;&gt;  DIST:1.638300000e+04  SSE_RUB:0.000000000e+00 SSE_TOT:4.944833778e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 4.945e+14   UNWHT_SSE: 4.945e+14   WHT_SD: 2.486e+05   UNWHT_SD: 2.486e+05
===============================
CALC JAC    15&gt;&gt;&gt;  DIST:1.638300000e+04
CALC FUN    16&gt;&gt;&gt;  DIST:3.276700000e+04  SSE_RUB:0.000000000e+00 SSE_TOT:4.608065586e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 4.608e+14   UNWHT_SSE: 4.608e+14   WHT_SD: 2.400e+05   UNWHT_SD: 2.400e+05
===============================
CALC JAC    16&gt;&gt;&gt;  DIST:3.276700000e+04
CALC FUN    17&gt;&gt;&gt;  DIST:6.553500000e+04  SSE_RUB:0.000000000e+00 SSE_TOT:4.063378221e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 4.063e+14   UNWHT_SSE: 4.063e+14   WHT_SD: 2.254e+05   UNWHT_SD: 2.254e+05
===============================
CALC JAC    17&gt;&gt;&gt;  DIST:6.553500000e+04
CALC FUN    18&gt;&gt;&gt;  DIST:1.310710000e+05  SSE_RUB:0.000000000e+00 SSE_TOT:3.489399568e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 3.489e+14   UNWHT_SSE: 3.489e+14   WHT_SD: 2.088e+05   UNWHT_SD: 2.088e+05
===============================
CALC JAC    18&gt;&gt;&gt;  DIST:1.310710000e+05
CALC FUN    19&gt;&gt;&gt;  DIST:1.530418700e+05  SSE_RUB:0.000000000e+00 SSE_TOT:3.450782038e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 3.451e+14   UNWHT_SSE: 3.451e+14   WHT_SD: 2.077e+05   UNWHT_SD: 2.077e+05
===============================
CALC JAC    19&gt;&gt;&gt;  DIST:1.530418700e+05
CALC FUN    20&gt;&gt;&gt;  DIST:1.530418700e+05  SSE_RUB:0.000000000e+00 SSE_TOT:3.450782038e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 3.451e+14   UNWHT_SSE: 3.451e+14   WHT_SD: 2.077e+05   UNWHT_SD: 2.077e+05
===============================
END FIT
9.095534 seconds elapsed, 20 func evals, 19 jac evals
</font>
</pre>

<b>As it can be seen from the fit, the result is not good since the generated data is not described by a constaint value. To improve the accuracy of the fit, let's change the code contained in the fitmodel.py by adding the three-dimensional polynomial:</b><br>

<pre><font color="#5FD7FF">import</font> os
<font color="#5FD7FF">from</font> functools <font color="#5FD7FF">import</font> <font color="#34E2E2"><b>reduce</b></font>

<font color="#5FD7FF">from</font> pysurfacefit.models.sympy <font color="#5FD7FF">import</font> ModelSympy
<font color="#5FD7FF">from</font> pysurfacefit.fitpars <font color="#5FD7FF">import</font> Par, Parameters

<font color="#5FD7FF">from</font> pysurfacefit.models.library <font color="#5FD7FF">import</font> generate_powers_layer_3d, poly3d

<font color="#FCE94F">class</font> <font color="#34E2E2"><b>fitmodel</b></font>(ModelSympy):
    <font color="#AD7FA8">&quot;&quot;&quot;</font>
<font color="#AD7FA8">    Test model for fitting data with just one parameter.</font>
<font color="#AD7FA8">    &quot;&quot;&quot;</font>
    <font color="#FCE94F">def</font> <font color="#34E2E2"><b>__init__</b></font>(self,calc_switch=<font color="#AD7FA8">&apos;numbified&apos;</font>):
        self.__check_symbolic__ = <font color="#34E2E2"><b>False</b></font>
        self.__calc_switch__ = calc_switch <font color="#34E2E2"># symbolic, lambdified, numbified</font>
        self.__components__ = {}

        <font color="#34E2E2"># Initialize empty parameters.</font>
        self.__params__ = Parameters()

        <font color="#34E2E2"># Add polynomial parameters</font>
        self.__powers__ = <font color="#34E2E2"><b>reduce</b></font>(<font color="#FCE94F">lambda</font> a,b:a+b,
            [generate_powers_layer_3d(n) <font color="#FCE94F">for</font> n <font color="#FCE94F">in</font> <font color="#34E2E2"><b>range</b></font>(<font color="#AD7FA8">1</font>,<font color="#AD7FA8">4</font>)]
        )

        <font color="#FCE94F">for</font> i,j,k <font color="#FCE94F">in</font> self.__powers__:
            self.__params__.append(group=<font color="#AD7FA8">&apos;poly&apos;</font>, pars=[
                Par(name=<font color="#AD7FA8">&apos;poly_%d_%d_%d&apos;</font>%(i,j,k), value=<font color="#AD7FA8">0.0</font>, flag=<font color="#34E2E2"><b>True</b></font>),
            ])

        <font color="#34E2E2"># Add constant parameter.</font>
        self.__params__.append(group=<font color="#AD7FA8">&apos;constant&apos;</font>, pars=[
            Par(name=<font color="#AD7FA8">&apos;constant&apos;</font>, value=<font color="#AD7FA8">0.0</font>, flag=<font color="#34E2E2"><b>True</b></font>),
        ])

        <font color="#FCE94F">def</font> <font color="#34E2E2"><b>__units__</b></font>(self):
        <font color="#FCE94F">return</font> {<font color="#AD7FA8">&quot;input&quot;</font>: {}, <font color="#AD7FA8">&quot;output&quot;</font>: <font color="#AD7FA8">&quot;None&quot;</font>}

    <font color="#FCE94F">def</font> <font color="#34E2E2"><b>__func__</b></font>(self,params,x,y,z):

        <font color="#34E2E2"># polynomial</font>
        p = params.get_values(group=<font color="#AD7FA8">&apos;poly&apos;</font>)
        poly = poly3d(p,x,y,z,self.__powers__)

        <font color="#34E2E2"># constant</font>
        constant = params[<font color="#AD7FA8">&apos;constant&apos;</font>].get_value()

        <font color="#34E2E2"># final result</font>
        res = poly + constant

        <font color="#34E2E2"># components</font>
        self.__components__[<font color="#AD7FA8">&apos;constant&apos;</font>] = constant

        <font color="#FCE94F">return</font> res

model = fitmodel()

<font color="#FCE94F">if</font> os.path.exists(<font color="#AD7FA8">&apos;fitmodel.csv&apos;</font>):
    model.load_params(<font color="#AD7FA8">&apos;fitmodel.csv&apos;</font>)
<font color="#FCE94F">else</font>:
    model.save_params(<font color="#AD7FA8">&apos;fitmodel.csv&apos;</font>)
</pre>

<b>Let's run the fit again with more advanced model:</b>

<pre>pysurfacefit --config config.ini --fit</pre><br>

<b>This will start the fitting process and will print the intermediate fit statistics to the STDOUT:</b><br>

<pre>
<font size="-2">
Treating sampledata fitgroup as FitPoints
BEGIN FIT
USING SCIPY.OPTIMIZE.LEAST_SQUARES: METHOD= trf
METHOD OPTIONS: {&apos;max_nfev&apos;: 200, &apos;bounds&apos;: [(-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf), (inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf)], &apos;jac&apos;: &lt;function Fitter.fit_least_squares.&lt;locals&gt;.&lt;lambda&gt; at 0x7f5def410560&gt;}

==========================================
&lt;&lt;&lt;&lt; calling __sympy_initialize_func__ &gt;&gt;&gt;&gt;
==========================================
Progress:
     - creating sympy objects for inputs
     - creating sympy objects for parameters
     - get the Sympy expression by calling function with the Sympy objects
     - create lambdified Python function from sympy expression
     - create compiled (numbified) code from the lambdified function
==========================================

CALC FUN     1&gt;&gt;&gt;  DIST:0.000000000e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.324527156e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.325e+14   UNWHT_SSE: 5.325e+14   WHT_SD: 2.580e+05   UNWHT_SD: 2.580e+05
===============================

==========================================
&lt;&lt;&lt;&lt; calling __sympy_initialize_jac__ &gt;&gt;&gt;&gt;
==========================================
Progress:
     - get the Sympy expression by calling function with the Sympy objects
     - create lambdified Python function from sympy expression
     - create compiled (numbified) code from the lambdified function
==========================================

CALC JAC     1&gt;&gt;&gt;  DIST:0.000000000e+00
CALC FUN     2&gt;&gt;&gt;  DIST:1.000000000e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.223291702e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.223e+14   UNWHT_SSE: 5.223e+14   WHT_SD: 2.555e+05   UNWHT_SD: 2.555e+05
===============================
CALC JAC     2&gt;&gt;&gt;  DIST:1.000000000e+00
CALC FUN     3&gt;&gt;&gt;  DIST:2.999995352e+00  SSE_RUB:0.000000000e+00 SSE_TOT:5.024183860e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.024e+14   UNWHT_SSE: 5.024e+14   WHT_SD: 2.506e+05   UNWHT_SD: 2.506e+05
===============================
CALC JAC     3&gt;&gt;&gt;  DIST:2.999995352e+00
CALC FUN     4&gt;&gt;&gt;  DIST:6.999923483e+00  SSE_RUB:0.000000000e+00 SSE_TOT:4.639397348e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 4.639e+14   UNWHT_SSE: 4.639e+14   WHT_SD: 2.408e+05   UNWHT_SD: 2.408e+05
===============================
CALC JAC     4&gt;&gt;&gt;  DIST:6.999923483e+00
CALC FUN     5&gt;&gt;&gt;  DIST:1.499903800e+01  SSE_RUB:0.000000000e+00 SSE_TOT:3.923332008e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 3.923e+14   UNWHT_SSE: 3.923e+14   WHT_SD: 2.215e+05   UNWHT_SD: 2.215e+05
===============================
CALC JAC     5&gt;&gt;&gt;  DIST:1.499903800e+01
CALC FUN     6&gt;&gt;&gt;  DIST:3.098657439e+01  SSE_RUB:0.000000000e+00 SSE_TOT:2.703038473e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 2.703e+14   UNWHT_SSE: 2.703e+14   WHT_SD: 1.838e+05   UNWHT_SD: 1.838e+05
===============================
CALC JAC     6&gt;&gt;&gt;  DIST:3.098657439e+01
CALC FUN     7&gt;&gt;&gt;  DIST:6.266232927e+01  SSE_RUB:0.000000000e+00 SSE_TOT:1.075696361e+14     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 1.076e+14   UNWHT_SSE: 1.076e+14   WHT_SD: 1.160e+05   UNWHT_SD: 1.160e+05
===============================
CALC JAC     7&gt;&gt;&gt;  DIST:6.266232927e+01
CALC FUN     8&gt;&gt;&gt;  DIST:1.160883053e+02  SSE_RUB:0.000000000e+00 SSE_TOT:1.093359850e+13     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 1.093e+13   UNWHT_SSE: 1.093e+13   WHT_SD: 3.697e+04   UNWHT_SD: 3.697e+04
===============================
CALC JAC     8&gt;&gt;&gt;  DIST:1.160883053e+02
CALC FUN     9&gt;&gt;&gt;  DIST:1.900605403e+02  SSE_RUB:0.000000000e+00 SSE_TOT:2.602526806e+12     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 2.603e+12   UNWHT_SSE: 2.603e+12   WHT_SD: 1.804e+04   UNWHT_SD: 1.804e+04
===============================
CALC JAC     9&gt;&gt;&gt;  DIST:1.900605403e+02
CALC FUN    10&gt;&gt;&gt;  DIST:3.975523481e+02  SSE_RUB:0.000000000e+00 SSE_TOT:1.507914353e+12     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 1.508e+12   UNWHT_SSE: 1.508e+12   WHT_SD: 1.373e+04   UNWHT_SD: 1.373e+04
===============================
CALC JAC    10&gt;&gt;&gt;  DIST:3.975523481e+02
CALC FUN    11&gt;&gt;&gt;  DIST:8.545961426e+02  SSE_RUB:0.000000000e+00 SSE_TOT:5.250210120e+11     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 5.250e+11   UNWHT_SSE: 5.250e+11   WHT_SD: 8.101e+03   UNWHT_SD: 8.101e+03
===============================
CALC JAC    11&gt;&gt;&gt;  DIST:8.545961426e+02
CALC FUN    12&gt;&gt;&gt;  DIST:1.639588009e+03  SSE_RUB:0.000000000e+00 SSE_TOT:1.170270161e+11     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 1.170e+11   UNWHT_SSE: 1.170e+11   WHT_SD: 3.825e+03   UNWHT_SD: 3.825e+03
===============================
CALC JAC    12&gt;&gt;&gt;  DIST:1.639588009e+03
CALC FUN    13&gt;&gt;&gt;  DIST:3.039288277e+03  SSE_RUB:0.000000000e+00 SSE_TOT:3.211110377e+10     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 3.211e+10   UNWHT_SSE: 3.211e+10   WHT_SD: 2.003e+03   UNWHT_SD: 2.003e+03
===============================
CALC JAC    13&gt;&gt;&gt;  DIST:3.039288277e+03
CALC FUN    14&gt;&gt;&gt;  DIST:6.369248010e+03  SSE_RUB:0.000000000e+00 SSE_TOT:2.288170913e+09     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 2.288e+09   UNWHT_SSE: 2.288e+09   WHT_SD: 5.348e+02   UNWHT_SD: 5.348e+02
===============================
CALC JAC    14&gt;&gt;&gt;  DIST:6.369248010e+03
CALC FUN    15&gt;&gt;&gt;  DIST:7.993910245e+03  SSE_RUB:0.000000000e+00 SSE_TOT:9.174559605e-16     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 9.175e-16   UNWHT_SSE: 9.175e-16   WHT_SD: 3.386e-10   UNWHT_SD: 3.386e-10
===============================
CALC JAC    15&gt;&gt;&gt;  DIST:7.993910245e+03
CALC FUN    16&gt;&gt;&gt;  DIST:7.993910245e+03  SSE_RUB:0.000000000e+00 SSE_TOT:1.549934571e-17     ==&gt; WEIGHTED_FIT
=====data group statistics=====
                    sampledata:   N: 8000   MIN_WHT: 1.0e+00   MAX_WHT: 1.0e+00   WHT_SSE: 1.550e-17   UNWHT_SSE: 1.550e-17   WHT_SD: 4.402e-11   UNWHT_SD: 4.402e-11
===============================
CALC JAC    16&gt;&gt;&gt;  DIST:7.993910245e+03
END FIT
11.802427 seconds elapsed, 16 func evals, 16 jac evals
</font>
</pre>


```python
%cd test

from fitmodel import model

model
```

    /home/roman/work/python/PyDev/PySurfaceFit/git-repo/pysurfacefit-master/showcase/test
    jeanny, Ver.3.0





      #  group     names                values  bounds         weights    flags
    ---  --------  ----------  ---------------  -----------  ---------  -------
      0  poly      poly_0_0_1   7500            (-inf, inf)          0        1
      1  poly      poly_0_1_0     -2.92517e-11  (-inf, inf)          0        1
      2  poly      poly_1_0_0      1.48667e-12  (-inf, inf)          0        1
      3  poly      poly_0_0_2  -2750            (-inf, inf)          0        1
      4  poly      poly_0_1_1     -3.5252e-13   (-inf, inf)          0        1
      5  poly      poly_0_2_0     -3.17951e-12  (-inf, inf)          0        1
      6  poly      poly_1_0_1     -1.55127e-13  (-inf, inf)          0        1
      7  poly      poly_1_1_0      2.3229e-13   (-inf, inf)          0        1
      8  poly      poly_2_0_0      1            (-inf, inf)          0        1
      9  poly      poly_0_0_3    300            (-inf, inf)          0        1
     10  poly      poly_0_1_2      7.33912e-15  (-inf, inf)          0        1
     11  poly      poly_0_2_1     -1.30851e-14  (-inf, inf)          0        1
     12  poly      poly_0_3_0      0.1          (-inf, inf)          0        1
     13  poly      poly_1_0_2      3.78018e-14  (-inf, inf)          0        1
     14  poly      poly_1_1_1      1.25425e-14  (-inf, inf)          0        1
     15  poly      poly_1_2_0      1.44824e-14  (-inf, inf)          0        1
     16  poly      poly_2_0_1      4.30556e-14  (-inf, inf)          0        1
     17  poly      poly_2_1_0     -1.16635e-14  (-inf, inf)          0        1
     18  poly      poly_3_0_0     -1.22386e-14  (-inf, inf)          0        1
     19  constant  constant       10            (-inf, inf)          0        1



<b>The code above lists the fitted parameters for the updated model.</b>


```python

```

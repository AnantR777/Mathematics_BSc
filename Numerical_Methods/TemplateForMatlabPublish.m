%% David Hewett
% MATH0033 Numerical Methods Computational homework 0
%
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)
% 
% This document is a template for you to use when submitting your
% computational homework. Please follow the exact formatting
% (with my name above replaced by your own name and computational
% homework 0 replaced by computational homework 1,2 or 3 as
% appropriate. The section marks %% split the code up into
% sections, each of which can be run individually by moving the
% cursor inside the section and hitting Ctrl-Enter on the
% keyboard. I recommend you use a separate section for each
% sub-part of an exercise, as this means that when you are
% debugging you only need to run the section of code you are
% working on at that time, rather than the whole file. Also, it means that 
% the outputs and figures will appear in the right place in the published 
% pdf, rather than all being outputted together at the end of the file. 
% To write comments and provide answers to the exercises, use the "%"
% symbol at the start of the line, as in this section here.
% Please be sure to explain your results in words using plenty of
% these comments, rather than just outputting numbers and figures
% without explanation, which I might not be able to interpret.
% Use axis labels (commands xlabel and ylabel) and figure titles
% (command title) to make it clear what your figures relate to.
% When you have finished your work you can use the "publish"
% command (there is a "publish" button in the "publish" tab at
% the top of the Matlab window) to output the code, your
% comments, your results and your figures, to a single file. When
% editing you may find the default html format the most useful as
% it updates quickly, but when you come to submit, please change
% to pdf format, by clicking on the little drop-down menu below
% the "publish" button, selecting "edit publishing options" and
% changing the output file format to pdf. Note that, confusingly, the 
% resulting pdf file will appear in a folder entitled "html". 
%
% (Note that if you like using latex then you can output to latex and compile the
% .tex file into a nicer pdf with better quality formulas - see the 
% more detailed comments at the end of the current file for more information. 
% But I should emphasize that the use of latex is completely optional - 
% I don't expect it.)
%%
% (b)
% 
% Here are some examples you might find helpful if you haven't
% used Matlab before. Matlab's main data structures are arrays.
% Here is how to construct a two-by-two matrix A and a two-by-one
% vector b, and then invert the system Ax=b using a direct method
% (LU factorization). Note that the ";" at the end of a line
% suppresses the output to the command line (i.e. stops Matlab 
% printing A and b to the display and just stores them internally).
A=[1,2;3,4];
b=[1;2];
x=A\b
%%
% (c)
% 
% To define a function you can use the anonymous function
% constructor "@". E.g. for a function of one variable
% f(x)=x^2-1, set
f=@(x)x.^2-1;
%%
% The .^2 rather than just ^2 is important as it allows f to be
% evaluated entrywise on an array, rather than just on a scalar.
% Similarly, entrywise multiplication of arrays needs x.*y rather
% than just x*y, as * without a dot means matrix multiplication.
f3=f(3)
x=linspace(-1,1,100);
fx=f(x);
plot(x,fx)
xlabel('x')
ylabel('f(x)')
%%
% Plotting is simple in Matlab and you just have to practice to
% get used to all the options.
%%
% *Exercise 2*
%
% (a)
%
% Functions of more than one variable can be defined similarly,
% using @:
g=@(x,y)3*x+x.*y-sin(y.^4)
g([1;2],[5;6])
%% 
% Note that this definition assumes that x and y are arrays of
% the same dimensions, and outputs an array of the same
% dimensions.
%%
% *Help and other tips*
%
% The Matlab Help documentation is very good and the "doc"
% command is your friend. E.g. type "doc sin" for information
% about how to use the function sin. There is also lots of help
% available on the internet in forums and online tutorials. The
% demo files (demo1, ...) I use in lectures will be posted on the
% Moodle page so you can copy syntax from there too. But if in doubt,
% google it!
%
% Here are some particular tips that might help you out.
%
% If you want to call one of the function files I supplied on the
% Moodle page (e.g. bisection.m) then you need to put a copy of
% the function file in your current Matlab working directory, or
% use "addpath" to tell Matlab where to find the codes if they
% are stored in a different folder. (But the former is easier I
% think.)
%
% Don't save your scripts with filenames which include spaces. So
% use "DaveHewett_Homework1" rather than "Dave Hewett Homework
% 1". You'll find that the latter causes problems - in particular
% when you use "publish" it won't generate figures.
%
% Another common publishing error is when you try to publish to a pdf file
% that you already have open in Adobe Acrobat Reader. Close the pdf and try
% again!
%
% Using sections (%%) you can just execute a section using
% Ctrl-Enter rather than running the whole code. But note that sometimes
% you might get error messages if you are referring to variables defined in
% a previous section that you haven't run yet.
%
% You can run a whole mfile by just pressing F5.
% 
% Before you try to "publish", make sure your code runs properly
% without generating any errors. (To check this, press F5 while 
% in the Matlab editor.) If not, it will not publish correctly.
%
% format long means answers are displayed in full precision (16
% digits). format compact means there's no extra space between
% output lines.
%
% Please put semicolons at the ends of your lines of Matlab code, as a 
% general rule. This suppresses the output to the command line. 
% Of course the exception is when you actually need to output results, e.g.
% if you want to output the solution to the problem or show iterations etc). 
% As an example, typing "x=linspace(-1,1,101);" will assign to the variable
% x an array containing 101 numbers equally spaced between -1 and 1. But if
% you omit the semicolon and type "x=linspace(-1,1,101)" then the same
% allocation will happen, but Matlab will print a list of all 101 numbers
% in the command window (and in your published pdf), which you (and I) 
% really don't need to see!
% 
% "zeros" and "ones" are useful for creating arrays.
%
% If you want Matlab to stop executing a command (e.g. because it is stuck
% in a long loop with pauses and you don't want to click through all of
% them, or you forgot to include a stopping criterion in your algorithm!)
% then click inside the Command Window and press Ctrl-c.
%
% sup norm of an array: norm(x,Inf) or max(abs(x)) 
%
% 2-norm of an array: norm(x) or sqrt(sum(x.^2))
%
% In a figure window, when in zoom mode you can right-click and
% select "horizontal" or "vertical" zoom. These are useful for
% inspecting plots. You can also click and drag to select an area
% to zoom in on.
%
% Make sure you understand the difference between , and ; when
% constructing arrays. Also, ' means conjugate transpose. Use .'
% for transpose without conjugation (these two are only different
% for complex-valued arrays of course).
%
% For a two-dimensional array x (i.e. a matrix), x(:,i) gives
% the ith column of x and x(i,:) gives the ith row of x.
%
% If doing the same procedure for a number of different parameter
% values, e.g. N=10,20,40,80, then don't copy and paste the same
% code multiple times. Use a for loop as in the following
% example:
%%
%  Nvec=[10,20,40,80];       % Array of parameter values
%  result=zeros(size(Nvec)); % Zero array of the same size as Nvec 
%                            % ready to receive the results you compute
%  for i=1:length(Nvec)
%       N=Nvec(i); 
%       ...                  % Insert the code to do computation using N 
%       result(i)=...        % Store result of the computation using N
%  end 
%  plot(Nvec, result)        % plot the results
%% 
% *Using latex commands in MATLAB Publish* 
%%
% You don't need to use latex commands (or even know what latex is) to get 
% full marks in the assignments - it's fine to just write formulas like 
% sqrt(x0)+x0^2 within your comments. But if you like using latex and want
% to do so then here are some tips. 
%
% The simplest way to use latex with publish is to just include standard 
% simple latex commands in your comments like $sqrt(x_0)+x_0^2$. When you 
% publish to pdf, MATLAB renders the formulae as images and includes them 
% in the right place in your file.
%
% Unfortunately the images are low resolution so the output can be hard to 
% read. If you want to output better quality latex formulas, then try 
% publishing to the ``latex'' file type rather than ``pdf'', then compile 
% the resulting .tex file using pdflatex. (You shouldn't need to edit the 
% .tex file but of course if you want to change formatting etc then you can
% do so.)
%% 
% *``Matlab'' latex*
%
% To get basic inline expressions like $f(x)=0$ you can use single dollar 
% signs, and to get basic equations you can use double dollar signs:
%
% $$ f(x)= f(x_0) + f'(x_0)(x-x_0) + \frac{f''(\xi)}{2}(x-x_0)^2. $$ 
%%
% *Pure latex*
% 
% To use more advanced latex commands (like begin(equation) etc) you can 
% embed sections of "pure" latex as in the following example. Note that 
% the following code does not get rendered at all if you publish to pdf, 
% but if you publish to latex then it will appear as a latex equation.
% 
% <latex>
% \begin{equation}
% f(x)= f(x_0) + f'(x_0)(x-x_0) + \frac{f''(\xi)}{2}(x-x_0)^2.  
% \end{equation}
% </latex>
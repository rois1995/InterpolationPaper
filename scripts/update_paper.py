#!/usr/bin/env python3
from pathlib import Path
import re

paper = Path(__file__).resolve().parents[1] / "paper.tex"
text = paper.read_text(encoding="utf-8")

# Preamble: plotting and table packages.
text = text.replace(
    r"\usepackage{array}\n\usepackage{enumitem}",
    r"\usepackage{array}\n\usepackage[dvipsnames]{xcolor}\n\usepackage{graphicx}\n\usepackage{siunitx}\n\usepackage{pgfplots}\n\usepackage{pgfplotstable}\n\usepgfplotslibrary{groupplots}\n\pgfplotsset{compat=1.18}\n\usepackage{enumitem}",
)

# Reusable paths and plot styles.
anchor = r"\newcommand{\calK}{\mathcal{K}}"
plot_setup = r'''
\newcommand{\DataDir}{data}
\newcommand{\GeneratedDir}{generated}

\pgfplotsset{
    transfer convergence axis/.style={
        xmode=log,
        ymode=log,
        x dir=reverse,
        grid=both,
        minor grid style={gray!15},
        major grid style={gray!30},
        xlabel={Effective spacing $h$},
        tick label style={font=\scriptsize},
        label style={font=\small},
        title style={font=\small},
        legend style={font=\scriptsize,draw=none},
        line width=0.8pt,
        mark size=2.0pt,
    },
    transfer bar axis/.style={
        grid=major,
        tick label style={font=\scriptsize},
        label style={font=\small},
        legend style={font=\scriptsize,draw=none},
    },
    scheme NN1/.style={black,solid,mark=*},
    scheme NN2/.style={NavyBlue,solid,mark=square*},
    scheme NN3/.style={ForestGreen,solid,mark=triangle*},
    scheme BAR2/.style={BrickRed,dashed,mark=o},
    scheme BAR3/.style={Orange,dashed,mark=square},
    scheme GP2/.style={Purple,densely dotted,mark=diamond*},
    scheme GP3/.style={Magenta,densely dotted,mark=triangle},
}

\sisetup{
    scientific-notation=true,
    round-mode=places,
    round-precision=3,
}
'''
if plot_setup.strip() not in text:
    text = text.replace(anchor, anchor + "\n" + plot_setup)

# Clarify that the limiter scales the complete Taylor increment at third order.
old_limiter = r'''Unlimited derivatives are used in the smooth manufactured-field tests so that the nominal interpolation order can be measured directly.
For transonic or otherwise under-resolved solutions, however, extrapolating a polynomial from neighboring values may create new extrema, negative density, or negative pressure.
A limiter may then scale the reconstructed increment before it is used by the transfer operator.
Classical choices include the Barth--Jespersen limiter~\cite{BarthJespersen1989}, the differentiable Venkatakrishnan limiter~\cite{Venkatakrishnan1995}, and higher-order accuracy-preserving limiter functions~\cite{MichalakOllivierGooch2009,Nishikawa2022}.
The limiter is treated as a separate stabilization choice rather than as part of the formal order definition: in smooth regions it should remain inactive, whereas close to a shock it may deliberately reduce a nominally high-order reconstruction toward a lower-order bounded one.
Its type, parameters, activation frequency, and any additional positivity correction are therefore reported separately.'''
new_limiter = r'''Unlimited derivatives are used in the smooth manufactured-field tests so that the nominal interpolation order can be measured directly.
For transonic or otherwise under-resolved solutions, however, extrapolating a polynomial from neighboring values may create new extrema, negative density, or negative pressure.
A node-based limiter may then scale the reconstructed Taylor increment before it is used by NN2, NN3, or BAR3.
For a limiter factor $\Phi_i\in[0,1]$, the second-order reconstruction becomes
\begin{equation}
    \phi_i+\Phi_i g_i\cdot(x-x_i),
\end{equation}
and the third-order reconstruction uses the same factor for both derivative terms,
\begin{equation}
    \phi_i
    +\Phi_i g_i\cdot(x-x_i)
    +\frac{1}{2}\Phi_i(x-x_i)^TH_i(x-x_i).
    \label{eq:limited-quadratic-reconstruction}
\end{equation}
Thus $\Phi_i\rightarrow0$ collapses the entire reconstruction to the bounded nodal value rather than leaving an unrestricted Hessian term active near a discontinuity.
Classical choices include the Barth--Jespersen limiter~\cite{BarthJespersen1989}, the differentiable Venkatakrishnan limiter~\cite{Venkatakrishnan1995}, and higher-order accuracy-preserving limiter functions~\cite{MichalakOllivierGooch2009,Nishikawa2022}.
The limiter is treated as a separate stabilization choice rather than as part of the formal order definition: in smooth regions it should remain inactive, whereas close to a shock it may deliberately reduce a nominally high-order reconstruction toward a lower-order bounded one.
Its type, parameters, activation frequency, and any additional positivity correction are therefore reported separately.
Galerkin projection does not possess a per-node Taylor increment and is bounded by the projection-level treatment introduced below.'''
text = text.replace(old_limiter, new_limiter)

# Add the projection-level limiter before mixed-element geometry.
projection_limiter = r'''
\subsubsection{Projection-level boundedness treatment}
\label{sec:galerkin-limiter}

The unconstrained $L_2$ projection is a global best approximation, but it does not satisfy a pointwise maximum principle.
Near a discontinuity it can therefore exhibit Gibbs-type over- and undershoots even when the projection remains conservative in the consistent sense.
Because GP2 and GP3 are obtained from the coupled system in Eq.~\eqref{eq:galerkin-system}, the node-based Taylor limiter of Eq.~\eqref{eq:limited-quadratic-reconstruction} cannot be applied directly.
A second-stage limiter is instead applied to the projected target field.

For every target vertex $j$, a local admissible interval $[\ell_j,u_j]$ is formed from the source-vertex values belonging to source simplices that overlap the target elements incident to $j$.
Using only source vertices prevents a reconstructed quadratic edge value from enlarging the admissible interval through its own Gibbs oscillation.
The projected value is clipped to
\begin{equation}
    \ell_j-\tau r
    \leq \phi_{B,j}
    \leq u_j+\tau r,
    \qquad
    r=\phi_A^{\max}-\phi_A^{\min},
    \label{eq:galerkin-soft-bound}
\end{equation}
where $\tau$ is a dimensionless deadband.
The hard limiter uses $\tau=0$ and enforces strict boundedness.
A small positive value gives the soft limiter, allowing smooth under-resolved extrema to exceed the discrete source range by a controlled amount.
The integral changed by clipping is subsequently redistributed among unsaturated target vertices using the bounded repair procedure of Section~\ref{sec:global-correction}.
For GP2 this preserves the lumped and consistent integrals simultaneously when the $P^1$ basis weights coincide with the finite-volume weights.
For vertex-only GP3, the operation preserves the selected lumped target integral but does not remove the separate contribution carried by discarded $P^2$ edge degrees of freedom unless an additional lumped-conservation repair is requested.
Boundedness and both conservation measures are consequently reported independently.
Bounded projection and repair strategies of this type follow the constrained Galerkin-transfer literature~\cite{FarrellMaddison2011,ShashkovWendroff2004}.
'''
marker = r"\subsection{Simplicial representation of mixed-element meshes}"
if projection_limiter.strip() not in text:
    text = text.replace(marker, projection_limiter + "\n" + marker)

# Expand static verification with the two actual manufactured cases.
static_marker = r'''Before any flow calculation, every base and globally corrected variant is verified on constant, affine, quadratic, and smooth non-polynomial manufactured fields in both two and three dimensions where available.'''
static_replacement = r'''Before any flow calculation, every base and globally corrected variant is verified on constant, affine, quadratic, and smooth non-polynomial manufactured fields in both two and three dimensions where available.
The primary smooth two-dimensional test is the asymmetric field
\begin{equation}
    \phi(x,y)
    =
    \rho_\infty
    \left[
        2+0.3\exp\!\left(0.8\frac{x}{L_x}+0.5\frac{y}{L_y}\right)
        \sin\!\left(2\pi\frac{x}{L_x}\right)
        \sin\!\left(2\pi\frac{y}{L_y}\right)
    \right].
    \label{eq:asymmetric-manufactured-field}
\end{equation}
The exponential envelope prevents positive and negative oscillatory lobes from cancelling in the domain integral, making the conservation defect a meaningful diagnostic on symmetric as well as irregular grids.
The discontinuous test is a two-plateau circular inclusion with density values $1.2$ and $3.6$ separated by a curved, non-mesh-aligned interface.
This case is used to measure the collapse of smooth-field convergence order, Gibbs over- and undershoots, limiter performance, and the interaction between boundedness and conservation.'''
text = text.replace(static_marker, static_replacement)

# Add automated data processing subsection before aerodynamic cases.
auto_subsection = r'''
\subsection{Automated convergence analysis and LaTeX plotting}
\label{sec:automated-postprocessing}

All static-transfer figures and order tables are generated from the comma-separated files distributed with the manuscript.
The script \texttt{scripts/build\_plot\_data.py} reads the complete long-form files \texttt{Comparison\_Errors.csv} and \texttt{Comparison\_Conservation.csv}, selects the density results for the asymmetric case, and writes compact wide tables for \texttt{pgfplots}.
It also reads the discontinuous-case files, recomputes convergence orders, and generates shortened labels for the boundedness comparison.
No numerical value is transcribed manually into the LaTeX source.

For consecutive levels, the pairwise observed order is
\begin{equation}
    p_{k\rightarrow k+1}
    =
    \frac{\log(e_k/e_{k+1})}{\log(h_k/h_{k+1})}.
    \label{eq:pairwise-observed-order}
\end{equation}
The script additionally reports a fitted order $p_{\mathrm{fit}}$ obtained as the slope of a least-squares line through all available points in the $(\log h,\log e)$ plane.
The fitted value is useful for the discontinuous test, for which the small three-level sequence contains visible pre-asymptotic scatter even though the expected $L_2$ rate is one half.
Both the fitted order and the finest-level pairwise order are retained in the generated files.

The plots in Section~\ref{sec:results} read the generated CSV tables directly through \texttt{pgfplots}.
The supplied \texttt{Makefile} declares the raw CSV files as dependencies of the generated data and the PDF.
Consequently, replacing any input CSV and running \texttt{make} recomputes the order tables and recompiles every affected plot automatically.
'''
marker2 = r"\subsection{Aerodynamic test cases}"
if auto_subsection.strip() not in text:
    text = text.replace(marker2, auto_subsection + "\n" + marker2)

# Replace the first three Results subsections with data-backed results.
results_start = text.index(r"\subsection{Formal accuracy and polynomial reproduction}")
results_end = text.index(r"\subsection{Native mesh sensitivity}")
results_block = r'''
\subsection{Smooth asymmetric manufactured field}
\label{sec:results-asymmetric}

Figure~\ref{fig:asym-l2-convergence} reports the relative density $L_2$ error for the asymmetric manufactured field of Eq.~\eqref{eq:asymmetric-manufactured-field}.
The source mesh is transferred to an independently generated target mesh with twice the nominal resolution, and the refinement sequence uses target resolutions $N=64$, $128$, and $256$.
The effective spacing stored in the CSV file is used on the abscissa, so the mixed-grid order does not rely on an assumed factor-of-two ratio.

\begin{figure}[htbp]
    \centering
    \begin{tikzpicture}
        \begin{groupplot}[
            group style={group size=2 by 1,horizontal sep=1.2cm},
            transfer convergence axis,
            width=0.45\textwidth,
            height=0.34\textwidth,
        ]
        \nextgroupplot[
            title={Structured quadrilateral grid},
            ylabel={Relative $L_2$ error},
            legend to name=legend:asym-l2,
            legend columns=4,
        ]
        \addplot[scheme NN1] table[col sep=comma,x=h_Struct,y=Struct_NN1]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{NN1}
        \addplot[scheme NN2] table[col sep=comma,x=h_Struct,y=Struct_NN2]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{NN2}
        \addplot[scheme NN3] table[col sep=comma,x=h_Struct,y=Struct_NN3]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{NN3}
        \addplot[scheme BAR2] table[col sep=comma,x=h_Struct,y=Struct_BAR2]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{BAR2}
        \addplot[scheme BAR3] table[col sep=comma,x=h_Struct,y=Struct_BAR3]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{BAR3}
        \addplot[scheme GP2] table[col sep=comma,x=h_Struct,y=Struct_GP2]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{GP2}
        \addplot[scheme GP3] table[col sep=comma,x=h_Struct,y=Struct_GP3]{\GeneratedDir/asymmetric_convergence_L2.csv}; \addlegendentry{GP3}

        \nextgroupplot[title={Unstructured mixed grid}]
        \addplot[scheme NN1] table[col sep=comma,x=h_Mixed,y=Mixed_NN1]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \addplot[scheme NN2] table[col sep=comma,x=h_Mixed,y=Mixed_NN2]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \addplot[scheme NN3] table[col sep=comma,x=h_Mixed,y=Mixed_NN3]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \addplot[scheme BAR2] table[col sep=comma,x=h_Mixed,y=Mixed_BAR2]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \addplot[scheme BAR3] table[col sep=comma,x=h_Mixed,y=Mixed_BAR3]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \addplot[scheme GP2] table[col sep=comma,x=h_Mixed,y=Mixed_GP2]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \addplot[scheme GP3] table[col sep=comma,x=h_Mixed,y=Mixed_GP3]{\GeneratedDir/asymmetric_convergence_L2.csv};
        \end{groupplot}
    \end{tikzpicture}
    \vspace{0.5em}
    \ref{legend:asym-l2}
    \caption{Relative density $L_2$ error for the smooth asymmetric transfer. The curves are read directly from \texttt{generated/asymmetric\_convergence\_L2.csv}.}
    \label{fig:asym-l2-convergence}
\end{figure}

\begin{table}[htbp]
    \centering
    \caption{Observed density $L_2$ convergence orders for the asymmetric test. The tabulated value is the slope fitted through all three refinement levels; the generated CSV also retains the finest pairwise order.}
    \label{tab:asym-orders}
    \pgfplotstabletypeset[
        col sep=comma,
        columns={scheme,nominal,Struct_fit,Mixed_fit},
        columns/scheme/.style={string type,column name={Scheme}},
        columns/nominal/.style={column name={Nominal},fixed,precision=0},
        columns/Struct_fit/.style={column name={Structured},fixed,precision=2},
        columns/Mixed_fit/.style={column name={Mixed},fixed,precision=2},
        every head row/.style={before row=\toprule,after row=\midrule},
        every last row/.style={after row=\bottomrule},
        columns/scheme/.style={string type,column type=l,column name={Scheme}},
        columns/nominal/.style={column type=c,column name={Nominal},fixed,precision=0},
        columns/Struct_fit/.style={column type=c,column name={Structured},fixed,precision=2},
        columns/Mixed_fit/.style={column type=c,column name={Mixed},fixed,precision=2},
    ]{\GeneratedDir/asymmetric_orders.csv}
\end{table}

All variants recover their nominal smooth-field order on both grid families.
The fitted rates are approximately one, two, and three for NN1--NN3; BAR2 and GP2 are second order; and BAR3 and GP3 exceed order three slightly over the tested range before approaching rates near $3.2$ between the two finest levels.
At $N=256$ on the mixed grid, BAR3 and GP3 reach relative density errors of approximately $2.6\times10^{-7}$ and $3.4\times10^{-7}$, respectively, whereas NN3 gives $2.6\times10^{-6}$.
The result confirms that the degree-two WLS reconstruction and the simplicial GP3 source lifting retain third-order behavior on the irregular mixed grid.

\subsection{Lumped conservation on the asymmetric field}
\label{sec:results-asymmetric-conservation}

Figure~\ref{fig:asym-conservation} shows the absolute defect in the nodal finite-volume integral before any optional global repair.
The asymmetric envelope is essential here: unlike a symmetric product of sines, it does not allow the integral error to vanish through cancellation.
The figure reports the solver-lumped metric of Eq.~\eqref{eq:fv-discrete-integral}; it must not be confused with the consistent finite-element integral preserved by the complete Galerkin field.

\begin{figure}[htbp]
    \centering
    \begin{tikzpicture}
        \begin{groupplot}[
            group style={group size=2 by 1,horizontal sep=1.2cm},
            transfer convergence axis,
            width=0.45\textwidth,
            height=0.34\textwidth,
        ]
        \nextgroupplot[
            title={Structured quadrilateral grid},
            ylabel={Absolute lumped-integral defect},
            legend to name=legend:asym-cons,
            legend columns=4,
        ]
        \addplot[scheme NN1] table[col sep=comma,x=h_Struct,y=Struct_NN1]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{NN1}
        \addplot[scheme NN2] table[col sep=comma,x=h_Struct,y=Struct_NN2]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{NN2}
        \addplot[scheme NN3] table[col sep=comma,x=h_Struct,y=Struct_NN3]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{NN3}
        \addplot[scheme BAR2] table[col sep=comma,x=h_Struct,y=Struct_BAR2]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{BAR2}
        \addplot[scheme BAR3] table[col sep=comma,x=h_Struct,y=Struct_BAR3]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{BAR3}
        \addplot[scheme GP2] table[col sep=comma,x=h_Struct,y=Struct_GP2]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{GP2}
        \addplot[scheme GP3] table[col sep=comma,x=h_Struct,y=Struct_GP3]{\GeneratedDir/asymmetric_conservation_defect.csv}; \addlegendentry{GP3}

        \nextgroupplot[title={Unstructured mixed grid}]
        \addplot[scheme NN1] table[col sep=comma,x=h_Mixed,y=Mixed_NN1]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \addplot[scheme NN2] table[col sep=comma,x=h_Mixed,y=Mixed_NN2]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \addplot[scheme NN3] table[col sep=comma,x=h_Mixed,y=Mixed_NN3]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \addplot[scheme BAR2] table[col sep=comma,x=h_Mixed,y=Mixed_BAR2]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \addplot[scheme BAR3] table[col sep=comma,x=h_Mixed,y=Mixed_BAR3]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \addplot[scheme GP2] table[col sep=comma,x=h_Mixed,y=Mixed_GP2]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \addplot[scheme GP3] table[col sep=comma,x=h_Mixed,y=Mixed_GP3]{\GeneratedDir/asymmetric_conservation_defect.csv};
        \end{groupplot}
    \end{tikzpicture}
    \vspace{0.5em}
    \ref{legend:asym-cons}
    \caption{Absolute density defect in the solver-lumped integral for the smooth asymmetric field. The optional global correction is not applied in this figure.}
    \label{fig:asym-conservation}
\end{figure}

On the structured mesh, GP2 reaches $4.3\times10^{-12}$ at the finest level, consistent with agreement between the $P^1$ basis-function integrals and the finite-volume lumping up to the linear-solver tolerance.
The same coincidence does not hold on the mixed mesh, where the finest GP2 lumped defect is $2.5\times10^{-7}$ even though the complete Galerkin field preserves its consistent integral.
GP3 exhibits a lumped defect of order $10^{-6}$--$10^{-7}$ because the restart state retains only vertex values and discards the edge degrees of freedom.
The raw NN and BAR defects decrease under refinement but are not identically zero.
Their globally corrected variants are evaluated separately because the correction changes the scalar integral after the local interpolation has already been completed.

\subsection{Discontinuous transfer: reduced order, boundedness, and conservation}
\label{sec:results-discontinuous}

Figure~\ref{fig:disc-l2-convergence} shows the relative density $L_2$ error for the two-plateau circular discontinuity.
The nominal smooth-field orders collapse because the error is concentrated in an $O(h)$ band around the unresolved interface.
The fitted rates in Table~\ref{tab:disc-orders} are close to the expected $L_2$ rate of one half: approximately $0.51$ on the structured grid and $0.43$--$0.45$ on the mixed grid.
The higher-order schemes reduce the error constant but do not recover their smooth-field asymptotic order at the jump.

\begin{figure}[htbp]
    \centering
    \begin{tikzpicture}
        \begin{groupplot}[
            group style={group size=2 by 1,horizontal sep=1.2cm},
            transfer convergence axis,
            width=0.45\textwidth,
            height=0.34\textwidth,
        ]
        \nextgroupplot[
            title={Structured quadrilateral grid},
            ylabel={Relative $L_2$ error},
            legend to name=legend:disc-l2,
            legend columns=4,
        ]
        \addplot[scheme NN1] table[col sep=comma,x=h_Struct,y=Struct_NN1]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{NN1}
        \addplot[scheme NN2] table[col sep=comma,x=h_Struct,y=Struct_NN2]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{NN2}
        \addplot[scheme NN3] table[col sep=comma,x=h_Struct,y=Struct_NN3]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{NN3}
        \addplot[scheme BAR2] table[col sep=comma,x=h_Struct,y=Struct_BAR2]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{BAR2}
        \addplot[scheme BAR3] table[col sep=comma,x=h_Struct,y=Struct_BAR3]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{BAR3}
        \addplot[scheme GP2] table[col sep=comma,x=h_Struct,y=Struct_GP2]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{GP2}
        \addplot[scheme GP3] table[col sep=comma,x=h_Struct,y=Struct_GP3]{\GeneratedDir/discontinuous_convergence_L2.csv}; \addlegendentry{GP3}

        \nextgroupplot[title={Unstructured mixed grid}]
        \addplot[scheme NN1] table[col sep=comma,x=h_Mixed,y=Mixed_NN1]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \addplot[scheme NN2] table[col sep=comma,x=h_Mixed,y=Mixed_NN2]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \addplot[scheme NN3] table[col sep=comma,x=h_Mixed,y=Mixed_NN3]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \addplot[scheme BAR2] table[col sep=comma,x=h_Mixed,y=Mixed_BAR2]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \addplot[scheme BAR3] table[col sep=comma,x=h_Mixed,y=Mixed_BAR3]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \addplot[scheme GP2] table[col sep=comma,x=h_Mixed,y=Mixed_GP2]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \addplot[scheme GP3] table[col sep=comma,x=h_Mixed,y=Mixed_GP3]{\GeneratedDir/discontinuous_convergence_L2.csv};
        \end{groupplot}
    \end{tikzpicture}
    \vspace{0.5em}
    \ref{legend:disc-l2}
    \caption{Relative density $L_2$ error for transfer of the discontinuous circular inclusion.}
    \label{fig:disc-l2-convergence}
\end{figure}

\begin{table}[htbp]
    \centering
    \caption{Fitted density $L_2$ convergence orders for the discontinuous transfer. The expected rate is $1/2$ independently of the nominal smooth-field order.}
    \label{tab:disc-orders}
    \pgfplotstabletypeset[
        col sep=comma,
        columns={scheme,nominal,Struct_fit,Mixed_fit},
        columns/scheme/.style={string type,column type=l,column name={Scheme}},
        columns/nominal/.style={column type=c,column name={Expected},fixed,precision=2},
        columns/Struct_fit/.style={column type=c,column name={Structured},fixed,precision=2},
        columns/Mixed_fit/.style={column type=c,column name={Mixed},fixed,precision=2},
        every head row/.style={before row=\toprule,after row=\midrule},
        every last row/.style={after row=\bottomrule},
    ]{\GeneratedDir/discontinuous_orders.csv}
\end{table}

At the finest structured level, the relative $L_2$ error ranges from approximately $6.97\times10^{-2}$ for NN1 to $5.43\times10^{-2}$ for BAR3 and $5.53\times10^{-2}$ for GP3.
The corresponding mixed-grid values are $6.58\times10^{-2}$, $4.92\times10^{-2}$, and $5.02\times10^{-2}$.
Thus the third-order methods improve the constant multiplying the discontinuity error, but the regularity of the transferred field determines the convergence rate.

Figure~\ref{fig:disc-boundedness} separates accuracy from admissibility.
Unlimited NN2, NN3, and BAR3 produce overshoots of $0.800$, $0.426$, and $0.133$, respectively, relative to the upper density plateau, accompanied by undershoots of $0.400$, $0.223$, and $0.140$.
The Barth--Jespersen limiter removes these extrema exactly in the tested configuration.
Applying the global conservation repair after limiting retains zero over- and undershoot while reducing the lumped conservation defect to $O(10^{-16})$.
The hard GP2 limiter is likewise strictly bounded and remains conservative to the linear-solver tolerance, whereas vertex-only GP3 remains bounded but retains a lumped defect of approximately $3.0\times10^{-4}$ associated with the discarded quadratic edge degrees of freedom.

\begin{figure}[htbp]
    \centering
    \begin{tikzpicture}
        \begin{axis}[
            transfer bar axis,
            width=0.96\textwidth,
            height=0.28\textwidth,
            ybar,
            bar width=5pt,
            ymin=0,
            ylabel={Bound violation},
            xtick=data,
            xticklabels=\empty,
            legend columns=2,
            legend style={at={(0.5,1.03)},anchor=south},
        ]
        \addplot[fill=NavyBlue!65,draw=NavyBlue] table[col sep=comma,x=idx,y=overshoot]{\GeneratedDir/discontinuous_boundedness_plot.csv};
        \addlegendentry{Overshoot}
        \addplot[fill=BrickRed!65,draw=BrickRed] table[col sep=comma,x=idx,y=undershoot]{\GeneratedDir/discontinuous_boundedness_plot.csv};
        \addlegendentry{Undershoot}
        \end{axis}
    \end{tikzpicture}

    \vspace{0.8em}
    \begin{tikzpicture}
        \begin{axis}[
            transfer bar axis,
            width=0.96\textwidth,
            height=0.30\textwidth,
            ybar,
            ymode=log,
            ymin=1e-17,
            ymax=2e-2,
            bar width=7pt,
            ylabel={Absolute conservation defect},
            xlabel={Configuration},
            xtick=data,
            xticklabels from table={\GeneratedDir/discontinuous_boundedness_plot.csv}{plot_label},
            x tick label style={rotate=55,anchor=east,font=\scriptsize},
        ]
        \addplot[fill=ForestGreen!65,draw=ForestGreen] table[col sep=comma,x=idx,y=abs_cons_defect]{\GeneratedDir/discontinuous_boundedness_plot.csv};
        \end{axis}
    \end{tikzpicture}
    \caption{Boundedness and lumped conservation for the discontinuous transfer. U denotes unlimited reconstruction, BJ the Barth--Jespersen limiter, GC the bounded global conservation correction, and S/H the soft/hard Galerkin projection limiter.}
    \label{fig:disc-boundedness}
\end{figure}

The refinement behavior of the uncorrected conservation defect is shown in Fig.~\ref{fig:disc-conservation}.
The structured GP2 curve remains at the numerical floor, while the mixed-grid lumped metric does not because the finite-volume weights and the consistent $P^1$ projection integral are not identical for that discretization.
This distinction confirms that conservation must always be stated together with the discrete integral being measured.

\begin{figure}[htbp]
    \centering
    \begin{tikzpicture}
        \begin{groupplot}[
            group style={group size=2 by 1,horizontal sep=1.2cm},
            transfer convergence axis,
            width=0.45\textwidth,
            height=0.34\textwidth,
        ]
        \nextgroupplot[
            title={Structured quadrilateral grid},
            ylabel={Absolute lumped-integral defect},
            legend to name=legend:disc-cons,
            legend columns=4,
        ]
        \addplot[scheme NN1] table[col sep=comma,x=h_Struct,y=Struct_NN1]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{NN1}
        \addplot[scheme NN2] table[col sep=comma,x=h_Struct,y=Struct_NN2]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{NN2}
        \addplot[scheme NN3] table[col sep=comma,x=h_Struct,y=Struct_NN3]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{NN3}
        \addplot[scheme BAR2] table[col sep=comma,x=h_Struct,y=Struct_BAR2]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{BAR2}
        \addplot[scheme BAR3] table[col sep=comma,x=h_Struct,y=Struct_BAR3]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{BAR3}
        \addplot[scheme GP2] table[col sep=comma,x=h_Struct,y=Struct_GP2]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{GP2}
        \addplot[scheme GP3] table[col sep=comma,x=h_Struct,y=Struct_GP3]{\GeneratedDir/discontinuous_conservation_defect.csv}; \addlegendentry{GP3}

        \nextgroupplot[title={Unstructured mixed grid}]
        \addplot[scheme NN1] table[col sep=comma,x=h_Mixed,y=Mixed_NN1]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \addplot[scheme NN2] table[col sep=comma,x=h_Mixed,y=Mixed_NN2]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \addplot[scheme NN3] table[col sep=comma,x=h_Mixed,y=Mixed_NN3]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \addplot[scheme BAR2] table[col sep=comma,x=h_Mixed,y=Mixed_BAR2]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \addplot[scheme BAR3] table[col sep=comma,x=h_Mixed,y=Mixed_BAR3]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \addplot[scheme GP2] table[col sep=comma,x=h_Mixed,y=Mixed_GP2]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \addplot[scheme GP3] table[col sep=comma,x=h_Mixed,y=Mixed_GP3]{\GeneratedDir/discontinuous_conservation_defect.csv};
        \end{groupplot}
    \end{tikzpicture}
    \vspace{0.5em}
    \ref{legend:disc-cons}
    \caption{Absolute solver-lumped conservation defect for the discontinuous transfer before optional global correction.}
    \label{fig:disc-conservation}
\end{figure}

'''
text = text[:results_start] + results_block + text[results_end:]

# Add the static findings to the Conclusions before the placeholder paragraph.
conclusion_anchor = r'''The same framework explains why a spatially remote transfer error may become force-visible only after acoustic or convective propagation.

\textit{[Replace this paragraph with quantitative conclusions after the numerical campaign.'''
conclusion_replacement = r'''The same framework explains why a spatially remote transfer error may become force-visible only after acoustic or convective propagation.

The completed static tests show that NN1--NN3 attain first-, second-, and third-order accuracy on both structured and mixed grids, while BAR2/GP2 and BAR3/GP3 attain their nominal second- and third-order behavior.
At a discontinuity, all families collapse to an $L_2$ rate close to one half, although higher-order transfer reduces the error constant.
Unlimited derivative-based point transfer is not bounded; the Barth--Jespersen limiter removes the measured extrema, and the subsequent global repair restores the lumped integral to machine precision without reintroducing them.
Hard-limited GP2 is simultaneously bounded and conservative in the tested configuration, whereas vertex-only GP3 retains the expected distinction between consistent projection conservation and the solver-lumped integral.

\textit{[Replace this paragraph with quantitative conclusions after the unsteady aerodynamic campaign.'''
text = text.replace(conclusion_anchor, conclusion_replacement)

paper.write_text(text, encoding="utf-8")
print(f"Updated {paper}")

(linear_uncertainty_propagation_framework)=
## Linear uncertainty propagation framework

The purpose of this internal document is to provide the theoretical foundations and the practical framework required to undertake an uncertainty analysis on a set of experimental measurements. Although uncertainty quantification is an essential part of any sound experimental methodology, it is often treated as a secondary task in day-to-day laboratory practice, where the main attention is usually devoted to instrumentation, acquisition, and post-processing. As a result, many measurement campaigns yield mean values, trends, and performance indicators without a sufficiently explicit assessment of how reliable those reported quantities are. The aim of the present document is to address that gap by presenting a structured route from the recorded signals to the corresponding uncertainty estimates.

The document is intended to be both instructive and operational. On the one hand, it develops the statistical and metrological ideas required to construct a rigorous uncertainty budget, clarifying the meaning of the quantities involved and the assumptions under which they are valid. On the other hand, it is meant to provide a systematic framework that can later be mapped onto a particular experimental set-up in a consistent and reproducible way. Although every measurement campaign is inherently case-specific, the logic connecting direct measurements, derived quantities, and propagated uncertainties remains general. For that reason, the present chapter first establishes the generic framework and only afterwards turns to the particular case of the wind-tunnel measurements.

The analysis is based on the linear uncertainty propagation framework introduced by Kline and McClintock {cite}`1953Kline`. In essence, that framework expresses the uncertainty of a derived quantity in terms of the uncertainties of the variables from which it is computed, weighted by the corresponding sensitivities of the functional relation. In practice, however, the successful application of that framework requires more than simply evaluating partial derivatives. One must first establish what the measured quantities represent, how their representative values are extracted from time records, how the random uncertainty of those representative values is estimated, and how non-random contributions are incorporated into the final uncertainty budget. The discussion therefore begins with the statistical treatment of steady measurement records.

---

(statistical_treatment_of_steady_measurement_records)=
## Statistical treatment of steady measurement records

In many experimental campaigns, the quantities of interest are associated with what may be termed a **static experimental configuration**. By this, one does not mean that the instantaneous signals are perfectly constant in time, but rather that the physical operating point of the system is intended to remain unchanged during the acquisition window. A given physical configuration is fixed, and the corresponding recorded signals fluctuate around a statistically steady mean value. In such cases, the purpose of the measurement is not to characterize a transient evolution, but to estimate one or several representative mean quantities associated with that nominal configuration.

Under such conditions, the raw outcome of the measurement system is typically a set of time series:

$$
x_1(t),x_2(t),\dots,x_n(t),
$$

each one corresponding to a different measured channel. The first task is therefore to determine which statistical quantities are extracted from those records, and what each of them means physically.

(the_mean)=
### The mean

For a given signal $x(t)$, sampled over a finite record of duration $T_{\mathrm{rec}}$ and containing $N$ samples, the first statistical quantity of interest is the sample mean:

$$
\overline{x}=\frac{1}{N}\sum_{i=1}^{N}x_i.
$$

This quantity is taken as the estimator of the mean value associated with the static operating condition. In most steady experimental measurements, this average is the quantity ultimately reported as the representative outcome of the acquisition.

:::{note}
In a wind-tunnel context, examples of such mean quantities include:
- the mean inlet velocity, $\overline{V}$,
- the mean ambient pressure, $\overline{P}$,
- the mean ambient temperature, $\overline{T}$,
- the mean load readings, $\overline{F}_{x,y,z}$,
- or the mean pressure deficit at a given wake-rake location, $\overline{p}_{i}$.
:::

(the_standard_deviation)=
### The standard deviation of the signal

A second quantity naturally obtained from the same record is the sample standard deviation:

$$
s_x=\sqrt{\frac{1}{N-1}\sum_{i=1}^{N}(x_i-\overline{x})^2}.
$$

This quantity measures the scatter of the **instantaneous** signal around its mean. It therefore describes the fluctuation level of the signal itself. If the record is very noisy or physically unsteady, $s_x$ will be large; if the signal fluctuates only weakly around the mean, $s_x$ will be small.

It is important to stress that $s_x$ is not, by itself, the uncertainty of the reported mean value. Rather, it is a measure of the spread of the recorded instantaneous values. This distinction is essential in steady experimental measurements, because the variability of the signal and the uncertainty of the estimated mean are not the same quantity.

(the_uncertainty_of_the_mean)=
### The uncertainty of the mean

Suppose that the objective is to report the mean value $\overline{x}$ associated with the nominal operating condition. The question is then no longer how strongly the instantaneous signal fluctuates, but how accurately the finite record estimates the true mean of the underlying process.

If the $N$ recorded samples were mutually independent, the standard uncertainty of the sample mean would be estimated as:

$$
\widehat{\sigma}_{\overline{x}}=\frac{s_x}{\sqrt{N}}.
$$

This is the classical **standard error of the mean**. It follows from the fact that the variance of the average of $N$ independent observations is the variance of the underlying variable divided by $N$.

At this point, the difference between the two quantities can already be stated clearly:

- the standard deviation $s_x$ characterizes the amplitude of the instantaneous fluctuations of the signal,
- whereas the standard uncertainty $\widehat{\sigma}_{\overline{x}}$ characterizes the precision with which the mean value has been estimated from the available record.

Thus, even before introducing any uncertainty-propagation framework, a steady measurement record already yields two distinct pieces of statistical information: the variability of the signal, and the uncertainty of the reported mean.

:::{important}
The standard deviation of the signal and the uncertainty of the mean are not the same quantity.

For a finite record of $N$ samples, the sample standard deviation is:

$$
s_x=\sqrt{\frac{1}{N-1}\sum_{i=1}^{N}(x_i-\overline{x})^2},
$$

and it estimates the fluctuation level of the recorded signal itself.

If the samples are independent, the estimated standard uncertainty of the mean is:

$$
\widehat{\sigma}_{\overline{x}}=\frac{s_x}{\sqrt{N}},
$$

which instead measures how accurately the finite record estimates the true mean of the underlying process.

Thus, $s_x$ and $\widehat{\sigma}_{\overline{x}}$ answer different questions:
- $s_x$ describes the spread of the instantaneous values,
- $\widehat{\sigma}_{\overline{x}}$ describes the precision of the reported average.
:::

(why_temporal_correlation_matters)=
### Why temporal correlation matters

The expression

$$
\widehat{\sigma}_{\overline{x}}=\frac{s_x}{\sqrt{N}}
$$

is valid only if the samples are statistically independent. In many real laboratory measurements, however, this assumption is not satisfied. Even under nominally steady operating conditions, time series may exhibit temporal correlations due to instrumentation dynamics, structural response, filtering, fluid-dynamic coherence, or slow facility-level variations. In such cases, consecutive samples do not necessarily provide fully independent pieces of information.

This is particularly relevant in wind-tunnel measurements. Pressure, velocity, and load records may retain memory over finite time intervals because turbulence structures persist, the wake exhibits coherent motion, pneumatic lines have response times, and measurement chains may behave as low-pass filters. If the sampling frequency is high, many successive samples may therefore be strongly correlated.

When this happens, the raw number of samples $N$ overestimates the actual amount of independent information contained in the record. The relevant quantity is then not $N$, but the **effective number of statistically independent samples**, denoted here by $N_{\mathrm{eff},x}$. The corresponding uncertainty of the mean must be written as:

$$
\widehat{\sigma}_{\overline{x},\mathrm{corr.}}=\frac{s_x}{\sqrt{N_{\mathrm{eff},x}}},
$$

with:

$$
N_{\mathrm{eff},x}<N
$$

whenever temporal correlation is present.

The natural next question is therefore how $N_{\mathrm{eff},x}$ is actually obtained from a measured time series. A convenient and physically meaningful route is provided by the normalized autocorrelation function.

---

(calculating_Neff)=
## Calculating $N_{\mathrm{eff},x}$

For a time signal $x(t)$, we first decompose the signal into its mean and fluctuating part:

$$
x(t)=\mu_x+x'(t),
$$

where:

$$
x'(t)=x(t)-\mu_x.
$$

The fluctuating component $x'(t)$ contains the departures of the signal from its mean value. Since the uncertainty of the mean is controlled by how those fluctuations persist in time, the relevant object is the covariance between the signal and a delayed copy of itself.

(on_the_normalized_autocorrelation)=
### The normalized autocorrelation

The autocovariance at lag $\tau$ is defined as:

$$
R_{xx}(\tau)=\mathbb E\left[(x(t)-\mu_x)(x(t+\tau)-\mu_x)\right]
=
\mathbb E\left[x'(t)x'(t+\tau)\right].
$$

Mathematically, $R_{xx}(\tau)$ is the covariance between the signal and a delayed copy of itself. Physically, it quantifies the memory of the fluctuation field. If a positive fluctuation at time $t$ tends to be followed by another positive fluctuation at time $t+\tau$, and a negative fluctuation tends to be followed by another negative fluctuation, then the product $x'(t)x'(t+\tau)$ is mostly positive and $R_{xx}(\tau)>0$. If the two fluctuations are unrelated, positive and negative products cancel out and $R_{xx}(\tau)\approx 0$. If the signal tends to swing to the opposite side after the delay $\tau$, then $R_{xx}(\tau)<0$.

At zero lag, the signal is compared with itself:

$$
R_{xx}(0)=\mathbb E\left[x'^2(t)\right]=\sigma_x^2.
$$

Therefore, the autocovariance starts from the variance of the signal and then decays as the signal loses memory of its previous state.

If the autocovariance is divided by the variance $R_{xx}(0)=\sigma_x^2$, the **normalized autocorrelation** results:

$$
\rho_{xx}(\tau)=\frac{R_{xx}(\tau)}{R_{xx}(0)}=\frac{R_{xx}(\tau)}{\sigma_x^2}.
$$

This normalization is what makes the function physically readable:

- $\rho_{xx}(0)=1$ always,
- $\rho_{xx}(\tau)\approx 1$ means the signal still remembers strongly where it was,
- $\rho_{xx}(\tau)\approx 0$ means the signal has essentially lost memory,
- $\rho_{xx}(\tau)<0$ means the signal tends to swing to the opposite side after that lag.

This is precisely why the autocorrelation is useful in the present context. The uncertainty of the mean depends not only on how large the fluctuations are, but also on how long they persist. The standard deviation $s_x$ measures the first of those aspects; the autocorrelation measures the second.

(the_integral_time_scale)=
### The integral time scale

The information contained in the normalized autocorrelation can be summarized through the **integral time scale**, defined in continuous form as:

$$
\tau_{\mathrm{int},x}=\int_0^\infty \rho_{xx}(\tau)\,\mathrm{d}\tau.
$$

This expression has a direct physical meaning. Since $\rho_{xx}(\tau)$ measures the fraction of fluctuation memory remaining after a delay $\tau$, the integral measures the total amount of memory contained in the signal. A signal whose autocorrelation decays rapidly has a small area and therefore a small integral time scale. A signal whose autocorrelation decays slowly has a larger area and therefore a larger integral time scale.

In practice, the integral cannot be extended to infinity. At large lags, the autocorrelation estimate becomes noisy because fewer overlapping data pairs are available. In addition, once the autocorrelation crosses zero, the dominant positive-memory contribution has usually ended. Therefore, one commonly truncates the integral at a finite lag:

$$
\tau_m=m\Delta t,
$$

often chosen as the last positive lag before the first non-positive value of $\rho_{xx}$.

The truncated integral is then:

$$
\tau_{\mathrm{int},x}\approx\int_0^{\tau_m}\rho_{xx}(\tau)\,\mathrm{d}\tau.
$$

For sampled data, the available autocorrelation values are not continuous. They are known at discrete lag times:

$$
\tau_k=k\Delta t,\qquad k=0,1,2,\dots,m.
$$

Therefore, the integral must be approximated numerically. Using the trapezoidal rule:

$$
\int_0^{\tau_m}\rho_{xx}(\tau)\,\mathrm{d}\tau
\approx
\sum_{k=0}^{m-1}\frac{\rho_{xx}(\tau_k)+\rho_{xx}(\tau_{k+1})}{2}\Delta t.
$$

Substituting $\tau_k=k\Delta t$ gives:

$$
\tau_{\mathrm{int},x}\approx\frac{\Delta t}{2}\sum_{k=0}^{m-1}\left[\rho_{xx}(k\Delta t)+\rho_{xx}((k+1)\Delta t)\right].
$$

Expanding the first few terms:

$$
\begin{aligned}
\tau_{\mathrm{int},x}
&\approx
\frac{\Delta t}{2}
\big[
\rho_{xx}(0)
+
\rho_{xx}(\Delta t)
\\
&\phantom{\approx \frac{\Delta t}{2}\big[}
+
\rho_{xx}(\Delta t)
+
\rho_{xx}(2\Delta t)
\\
&\phantom{\approx \frac{\Delta t}{2}\big[}
+
\rho_{xx}(2\Delta t)
+
\rho_{xx}(3\Delta t)
\\
&\phantom{\approx \frac{\Delta t}{2}\big[}
+
\dots
+
\rho_{xx}((m-1)\Delta t)
+
\rho_{xx}(m\Delta t)
\big].
\end{aligned}
$$

Collecting repeated terms:

$$
\tau_{\mathrm{int},x}\approx\Delta t\left[\frac{1}{2}\rho_{xx}(0)+\sum_{k=1}^{m-1}\rho_{xx}(k\Delta t)+\frac{1}{2}\rho_{xx}(m\Delta t)\right].
$$

The reason the zero-lag term receives a weight of $\frac{1}{2}$ is simple: it is one of the two endpoints of the integration interval. In the trapezoidal rule, the endpoints are counted with half weight, while the interior points are counted with full weight.

Since:

$$
\rho_{xx}(0)=1,
$$

the expression becomes:

$$
\tau_{\mathrm{int},x}\approx\Delta t\left[\frac{1}{2}+\sum_{k=1}^{m-1}\rho_{xx}(k\Delta t)+\frac{1}{2}\rho_{xx}(m\Delta t)\right].
$$

In the implementation used in the script, a closely related approximation is adopted:

$$
\tau_{\mathrm{int},x}\approx\Delta t\left[\frac{1}{2}+\sum_{k=1}^{m}\rho_{xx}(k\Delta t)\right].
$$

This expression keeps the half-weighted zero-lag contribution but gives full weight to the last included positive lag. Since $m$ is usually selected near the first zero crossing, $\rho_{xx}(m\Delta t)$ is typically already small. Therefore, the practical difference between the strict trapezoidal expression and the implemented expression is usually minor.

The quantity $\tau_{\mathrm{int},x}$ is therefore the area under the normalized autocorrelation curve while it remains meaningfully positive. Physically, it is the typical duration over which one fluctuation still influences subsequent measurements.

(from_integral_time_scale_to_Neff)=
### From the integral time scale to $N_{\mathrm{eff},x}$

The next step is to relate the integral time scale to the uncertainty of the mean.

For a continuous record of duration $T_{\mathrm{rec}}$, the time mean is:

$$
\overline{x}=\frac{1}{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}}x(t)\,\mathrm{d}t.
$$

Since only the fluctuating part contributes to the uncertainty of the mean, we can write:

$$
\overline{x}-\mu_x=\frac{1}{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}}x'(t)\,\mathrm{d}t.
$$

The variance of the mean is therefore:

$$
\operatorname{Var}(\overline{x})=
\mathbb E\left[
\left(
\frac{1}{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}}x'(t)\,\mathrm{d}t
\right)^2
\right].
$$

Taking the square inside the expectation gives:

$$
\operatorname{Var}(\overline{x})=
\frac{1}{T_{\mathrm{rec}}^2}
\mathbb E\left[
\int_0^{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}}
x'(t_1)x'(t_2)\,\mathrm{d}t_1\mathrm{d}t_2
\right].
$$

Interchanging expectation and integration:

$$
\operatorname{Var}(\overline{x})=
\frac{1}{T_{\mathrm{rec}}^2}
\int_0^{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}}
\mathbb E\left[x'(t_1)x'(t_2)\right]
\,\mathrm{d}t_1\mathrm{d}t_2.
$$

For a stationary signal, the covariance between two points depends only on their time separation, not on the absolute time. Therefore:

$$
\mathbb E\left[x'(t_1)x'(t_2)\right]=R_{xx}(t_2-t_1).
$$

Hence:

$$
\operatorname{Var}(\overline{x})=
\frac{1}{T_{\mathrm{rec}}^2}
\int_0^{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}}
R_{xx}(t_2-t_1)\,\mathrm{d}t_1\mathrm{d}t_2.
$$

This double integral accounts for the covariance between every pair of points in the record. To simplify it, introduce the lag variable:

$$
\tau=t_2-t_1.
$$

For a given lag $\tau$, the number of pairs of points separated by that lag is proportional to:

$$
T_{\mathrm{rec}}-|\tau|.
$$

Therefore, the double integral over the square domain $0\leq t_1,t_2\leq T_{\mathrm{rec}}$ can be rewritten as an integral over lags:

$$
\int_0^{T_{\mathrm{rec}}}\int_0^{T_{\mathrm{rec}}} R_{xx}(t_2-t_1)\,\mathrm{d}t_1\mathrm{d}t_2 = \int_{-T_{\mathrm{rec}}}^{T_{\mathrm{rec}}}
\left(T_{\mathrm{rec}}-|\tau|\right)R_{xx}(\tau)\,\mathrm{d}\tau.
$$

Thus:

$$
\operatorname{Var}(\overline{x})=
\frac{1}{T_{\mathrm{rec}}^2}
\int_{-T_{\mathrm{rec}}}^{T_{\mathrm{rec}}}
\left(T_{\mathrm{rec}}-|\tau|\right)R_{xx}(\tau)\,\mathrm{d}\tau.
$$

Because the autocovariance of a real-valued stationary signal is an even function,

$$
R_{xx}(-\tau)=R_{xx}(\tau),
$$

we can write:

$$
\operatorname{Var}(\overline{x})=
\frac{2}{T_{\mathrm{rec}}^2}
\int_0^{T_{\mathrm{rec}}}
\left(T_{\mathrm{rec}}-\tau\right)R_{xx}(\tau)\,\mathrm{d}\tau.
$$

Factoring out $T_{\mathrm{rec}}$:

$$
\operatorname{Var}(\overline{x})=
\frac{2}{T_{\mathrm{rec}}}
\int_0^{T_{\mathrm{rec}}}
\left(1-\frac{\tau}{T_{\mathrm{rec}}}\right)R_{xx}(\tau)\,\mathrm{d}\tau.
$$

For a sufficiently long record, the autocovariance decays to zero over a time much smaller than the total record duration. That is:

$$
\tau_{\mathrm{int},x}\ll T_{\mathrm{rec}}.
$$

Over the range where $R_{xx}(\tau)$ is meaningfully nonzero:

$$
1-\frac{\tau}{T_{\mathrm{rec}}}\approx 1.
$$

Also, because $R_{xx}(\tau)$ is negligible at large lags, the upper limit can be extended from $T_{\mathrm{rec}}$ to $\infty$:

$$
\operatorname{Var}(\overline{x})\approx\frac{2}{T_{\mathrm{rec}}}\int_0^\infty R_{xx}(\tau)\,\mathrm{d}\tau.
$$

Using:

$$
R_{xx}(\tau)=\sigma_x^2\rho_{xx}(\tau),
$$

we obtain:

$$
\operatorname{Var}(\overline{x})\approx
\frac{2\sigma_x^2}{T_{\mathrm{rec}}}
\int_0^\infty \rho_{xx}(\tau)\,\mathrm{d}\tau.
$$

Since:

$$
\tau_{\mathrm{int},x}=\int_0^\infty \rho_{xx}(\tau)\,\mathrm{d}\tau,
$$

then:

$$
\operatorname{Var}(\overline{x})
\approx
\frac{2\sigma_x^2\tau_{\mathrm{int},x}}{T_{\mathrm{rec}}}.
$$

For an equivalent set of $N_{\mathrm{eff},x}$ independent samples, the variance of the mean would be:

$$
\operatorname{Var}(\overline{x})=\frac{\sigma_x^2}{N_{\mathrm{eff},x}}.
$$

Equating both expressions gives:

$$
\frac{\sigma_x^2}{N_{\mathrm{eff},x}}=
\frac{2\sigma_x^2\tau_{\mathrm{int},x}}{T_{\mathrm{rec}}}.
$$

After cancelling $\sigma_x^2$:

$$
N_{\mathrm{eff},x}\approx\frac{T_{\mathrm{rec}}}{2\tau_{\mathrm{int},x}}.
$$

This is the key physical correction. If the signal loses memory very quickly, $\tau_{\mathrm{int},x}$ is small, $N_{\mathrm{eff},x}$ is large, and the mean is estimated very well. If the signal is slowly varying, $\tau_{\mathrm{int},x}$ is large, $N_{\mathrm{eff},x}$ is small, and the recorded signal contains only a limited number of independent realizations.

:::{important}
From the viewpoint of the uncertainty of the mean, a correlated record behaves as if it contained approximately $N_{\mathrm{eff},x}$ independent blocks of information, each with a characteristic duration of order $2\tau_{\mathrm{int},x}$.

Thus, the correlation-corrected standard uncertainty of the mean is:

$$
\widehat{\sigma}_{\overline{x},\mathrm{corr.}}=
\frac{s_x}{\sqrt{N_{\mathrm{eff},x}}},
$$

not simply $s_x/\sqrt{N}$.
:::

::::{important} **Worked example: computing $N_{\mathrm{eff},x}$ in a time-series**

A very concrete example helps. Suppose the sampling frequency of a signal is $1000\ \mathrm{Hz}$ and $T_{\mathrm{rec}}=60\ \mathrm{s}$. Then the raw number of samples is:

$$
N=60000.
$$

If the signal were white noise, the autocorrelation would drop to zero almost immediately, $\tau_{\mathrm{int},x}$ would be on the order of $\Delta t$, and $N_{\mathrm{eff},x}$ would be close to $N$. But if the signal remains correlated for, say, $0.05\ \mathrm{s}$, then:

$$
N_{\mathrm{eff},x}\approx\frac{60}{2\times 0.05}=600.
$$

So despite having $60000$ raw samples, the relevant information for estimating the mean is equivalent to about $600$ independent samples. That is a huge difference, and it is exactly why the normalized autocorrelation is needed.

The figure below shows the first 2 seconds of the mentioned raw time signal of 60 seconds sampled at a frequency of 1000 Hz.

:::{figure} Figs/01_raw_time_signal.svg
:name: fig-raw-time-signal
:width: 85%
:align: center
First 2 seconds of a raw time signal of 60 seconds sampled at a frequency of 1000 Hz.
:::

Upon such a figure, computing the normalized autocorrelation provides the shape shown in the illustration below. The shaded area corresponds to the integral time scale $\tau_{\mathrm{int},x}$.

:::{figure} Figs/02_normalized_autocorrelation.svg
:name: fig-normalized-autocorrelation
:width: 85%
:align: center
Normalized autocorrelation function of the above signal, with the shaded area showing the integral time scale ($\tau_{\mathrm{int},x}$) computation.
:::

Based on such a $\tau_{\mathrm{int},x}$ value, the original 2-seconds-long time series can be chunked into pieces of $2\tau_{\mathrm{int},x}$, namely the statistically independent data chunks. Such a division is shown in the figure below.

:::{figure} Figs/03_chunked_signal.svg
:name: fig-chunked-signal
:width: 85%
:align: center
The original 2-seconds-long time signal chunked into pieces of $2\tau_{\mathrm{int},x}$, namely the statistically independent data chunks.
:::

::::

---

(random_and_nonrandom_uncertainty_contributions)=
## Random and non-random uncertainty contributions

Once the uncertainty of the mean has been established for a steady record, the next step is to place that quantity inside a broader uncertainty budget. The uncertainty associated with a reported measurand generally contains two conceptually distinct parts: one associated with the random scatter of the observations, and another associated with non-random or bias-like effects that are not inferred from the scatter of the current record itself.

(typeA_contributions)=
### Type A contributions

Type A contributions are those evaluated from the statistical analysis of repeated observations. In the present framework, they are associated with the random part of the measurement process and are obtained from the recorded time series themselves. For a mean quantity extracted from a steady record, the Type A standard uncertainty is therefore linked to the uncertainty of the estimated mean, not simply to the standard deviation of the raw signal. Thus:

$$
u_A(\overline{x})=\frac{s_x}{\sqrt{N_{\mathrm{eff},x}}},
$$

where the independent-sample case is recovered as a special case when $N_{\mathrm{eff},x}=N$.

(typeB_contributions)=
### Type B contributions

Type B contributions, on the other hand, are evaluated from information other than the statistical scatter of the current record. They represent uncertainty components associated with calibration, finite resolution, specified accuracy limits, drift, hysteresis, linearity limits, reference standards, signal-conversion chains, or engineering assumptions embedded in the measurement model.

In practice, these contributions are often estimated from manufacturer documentation, calibration certificates, previous laboratory experience, reference-comparison tests, or bounded tolerances. They are not identified by the mere fact that they are “systematic” in a colloquial sense, but by the way in which they are evaluated.

(interpreting_calibration_sheets_and_probability_models)=
### Interpreting calibration sheets, specifications, and probability models

The passage from device specifications to uncertainty values is one of the most delicate parts of an experimental uncertainty analysis. A calibration certificate, a manufacturer datasheet, a stated tolerance, or a resolution limit does **not** by itself provide a standard uncertainty unless the information is interpreted through an appropriate probability model. In the language of Type B evaluation, this interpretation is based on scientific judgment using all relevant information available, including previous experience with the instrument, manufacturer specifications, calibration reports, reference data, and any auxiliary checks carried out in the laboratory.

At this stage of the uncertainty analysis, the role of the probability model is to convert the information provided by the instrument or by its documentation into a **standard uncertainty contribution**. In other words, one starts from an interval, tolerance, or bounded specification supplied externally, and translates it into the standard-uncertainty scale on which all individual contributions can later be combined. This is therefore a “backward” use of the probability model: one moves from an available interval to the corresponding standard deviation-like quantity.

Suppose that the available information on a given device is expressed through a symmetric bound:

$$
\pm a.
$$

The conversion of that bound into a standard uncertainty depends on the probability model assigned to the underlying error.

If the quoted value already corresponds to an expanded uncertainty associated with a known coverage factor $k$, then the corresponding standard uncertainty is:

$$
u=\frac{a}{k}.
$$

If, instead, the only information available is that the error lies somewhere between $-a$ and $+a$, with no reason to prefer central values over edge values, a **rectangular distribution** is typically adopted, yielding:

$$
u=\frac{a}{\sqrt{3}}.
$$

If values near zero are considered more plausible than values close to the bounds, a **triangular distribution** may be more appropriate, in which case:

$$
u=\frac{a}{\sqrt{6}}.
$$

These factors therefore serve to place all Type B contributions on a common standard-uncertainty basis. At this point, the purpose is not yet to report a 95\% interval for the final measurand, but to express each individual contribution in a form that is compatible with the uncertainty-propagation framework.

It is nevertheless useful to anticipate that the same probability models will later appear again, but in the opposite direction. Once the individual contributions have been converted into standard uncertainties and combined, the final result will be reported through a coverage interval, typically at the 95\% level. The same distributional assumptions then determine how one passes from a standard uncertainty back to an interval. In this sense, the uncertainty analysis uses probability models in two complementary ways:
- first, to move from specification intervals to standard uncertainties for the individual contributions;
- later, to move from the final combined standard uncertainty to the reported coverage interval.

This distinction is illustrated in the figure below, where the Gaussian, rectangular, and triangular distributions are shown together with their standard uncertainties and central 95\% coverage intervals. The figure makes visually clear that a factor such as $\sqrt{3}$ or $\sqrt{6}$ is not, in itself, a generic 95\% coverage factor. Rather, it is the factor that converts a bounded interval into the associated standard uncertainty under a given distributional assumption.

:::{figure} Figs/04_distributions_coverage_factors.svg
:name: fig-chunked-signal
:width: 85%
:align: center
Generic Gaussian (top), rectangular (middle) and triangular (bottom) distributions with their corresponding 95% coverage intervals.
:::

For a symmetric rectangular distribution with support $[-a,a]$, the standard uncertainty is $a/\sqrt{3}$, but the central 95\% interval is not obtained by multiplying $u$ by $\sqrt{3}$. Rather, since 95\% of the distribution occupies 95\% of the total width, the corresponding half-width is $0.95a$, and therefore:

$$
k_{95,\mathrm{rect.}}=\frac{0.95a}{a/\sqrt{3}}=0.95\sqrt{3}\approx 1.65.
$$

Similarly, for a symmetric triangular distribution with support $[-a,a]$, the standard uncertainty is $a/\sqrt{6}$, while the central 95\% interval yields:

$$
k_{95,\mathrm{tri.}}\approx 1.90.
$$

Thus, factors such as $\sqrt{3}$ or $\sqrt{6}$ are **conversion factors to standard uncertainty**, not generic 95\% coverage factors. This distinction should be made explicit whenever calibration sheets, tolerance bands, or probe resolutions are interpreted within the uncertainty budget.

This point is especially relevant when dealing with probe/device resolutions and with signal-conditioning chains. A finite resolution, a least displayed digit, or a bounded quantization step is typically treated as a bounded interval and therefore often modeled through a rectangular distribution. By contrast, a calibration certificate may already provide an expanded uncertainty associated with a stated coverage factor, and should then be converted back to a standard uncertainty through division by that factor. The same numerical bound may therefore lead to different standard uncertainties depending on what the bound actually represents. The role of the present framework is precisely to make those distinctions explicit and systematic, so that the assignment of Type B contributions is not left to ad hoc judgment.

---

(on_the_linear_propagation_framework)=
## Linear propagation of uncertainty

Once the individual uncertainty contributions associated with the input quantities have been established, they must be propagated through the measurement equations linking those inputs to the reported measurands. This is the point at which the Kline–McClintock framework enters explicitly.

(basic_and_derived_measurands)=
### Basic and derived measurands

A central distinction in this framework is the one between **basic measurands** and **derived measurands**.

A basic measurand is a quantity directly obtained from the measurement system, such as a pressure signal delivered by a pressure transducer, a temperature signal recorded by a thermometer, or a load signal read by a balance.

A derived measurand, by contrast, is obtained from a functional relationship involving one or more basic measurands, and possibly other already derived quantities. Density computed from pressure and temperature, dynamic pressure computed from velocity, Reynolds number computed from velocity, density, viscosity, and chord length, or aerodynamic coefficients computed from loads and reference scales are all examples of derived measurands.

The uncertainty analysis is therefore naturally hierarchical: one begins with the uncertainty characterization of the basic measurands and then propagates those contributions through the successive derived relations.

(the_uncertainty_propagation_expression)=
### The uncertainty propagation expression

Let a derived measurand $Y$ be defined as:

$$
Y=f(X_1,X_2,\dots,X_n),
$$

where the variables $X_1,\dots,X_n$ represent the basic input quantities. If the uncertainties affecting those input variables are sufficiently small for a first-order Taylor expansion to be adequate, then the deviation of $Y$ from its nominal value may be approximated as:

$$
\delta Y \approx
\frac{\partial f}{\partial X_1}\delta X_1+
\frac{\partial f}{\partial X_2}\delta X_2+
\dots+
\frac{\partial f}{\partial X_n}\delta X_n.
$$

The coefficients:

$$
\frac{\partial f}{\partial X_i}
$$

are the sensitivity coefficients of the derived quantity with respect to each input variable. They measure how strongly the result responds to a perturbation in each input, and they are the key bridge between the physical measurement equation and the uncertainty budget.

If the input quantities are treated as uncorrelated, the variance-like combination of those contributions leads to the familiar root-sum-square expression:

$$
u_Y^2 \approx
\left(\frac{\partial f}{\partial X_1}\right)^2u_{X_1}^2+
\left(\frac{\partial f}{\partial X_2}\right)^2u_{X_2}^2+
\dots+
\left(\frac{\partial f}{\partial X_n}\right)^2u_{X_n}^2,
$$

or, equivalently:

$$
u_Y \approx
\sqrt{
\sum_{i=1}^{n}
\left(\frac{\partial f}{\partial X_i}u_{X_i}\right)^2 }.
$$

If relevant correlations exist among the inputs, covariance terms must also be included. The key message remains the same: the uncertainty of a derived quantity is not assigned directly, but built systematically from the uncertainties of the quantities on which it depends.

---

(towards_expanded_uncertainty)=
## Reporting uncertainty intervals

At this stage, the uncertainty framework is naturally expressed in terms of **standard uncertainties**. This is the level at which Type A and Type B contributions are combined, and it is also the level at which the linear propagation law is formulated. Experimental results, however, are often reported not only through a best estimate and its associated standard uncertainty, but through an interval intended to contain the measurand with a prescribed level of coverage probability. In the terminology of the *Guide to the Expression of Uncertainty in Measurement* (GUM){cite}`2020GUM`, this is done through the **expanded uncertainty**.

In the original engineering spirit of Kline and McClintock{cite}`1952Kline`, input uncertainties were often stated on a common odds or confidence basis. In the present framework, however, all input contributions have been first converted to standard uncertainties and propagated on that common scale; a common 95% reporting interval is then recovered at the output level through the combined standard uncertainty, the effective degrees of freedom, and the corresponding distribution-based coverage factor.

(coverage_intervals)=
### Coverage intervals and expanded uncertainty

If the estimate of a measurand is denoted by $y$, and its combined standard uncertainty by $u_c(y)$, the expanded uncertainty associated with a target coverage probability $p$ is written as:

$$
U_p=k_p\,u_c(y),
$$

where $k_p$ is the coverage factor associated with the desired coverage level. The reported result is then expressed in the form:

$$
Y=y\pm U_p.
$$

Strictly speaking, GUM distinguishes between a **coverage interval** and the more specifically statistical notion of a **confidence interval**. In much experimental practice, however, the latter expression is used informally to refer to the former. In the present framework, the important point is that the standard uncertainty $u_c(y)$ quantifies uncertainty on a standard-deviation-like scale, whereas $U_p$ defines the interval actually reported around the estimate.

A particularly common choice in experimental work is the **95\% level**. When the distribution of the output quantity can be approximated as normal and the effective number of degrees of freedom is sufficiently large, one often uses:

$$
U_{95}\approx 2u_c(y).
$$

More precisely, for a normal distribution the central 95\% interval corresponds to:

$$
k_{95}=1.96,
$$

while the engineering approximation $k=2$ corresponds to a coverage probability slightly above 95\%. This approximation is widely used, but it should be understood as a large-sample simplification rather than as a universal rule.

(effective_degrees_of_freedom)=
### Effective degrees of freedom

The previous subsection dealt with the interpretation of externally provided intervals and bounds at the level of the **individual uncertainty contributions**. Once those contributions have been converted into standard uncertainties and propagated through the measurement model, the direction of reasoning is reversed. The task is no longer to infer a standard uncertainty from a given interval, but to determine which interval should be reported for the final measurand from the combined standard uncertainty $u_c(y)$.

In order to report a 95\% interval, however, one must also quantify how reliable that estimate of standard uncertainty is. This is done through the associated number of **degrees of freedom**.

For a standard uncertainty obtained from the standard deviation of the mean of $n$ independent repeated observations, the corresponding degrees of freedom are:

$$
\nu=n-1.
$$

When the number of degrees of freedom is finite, the appropriate 95\% coverage factor is taken from the Student $t$ distribution rather than replaced automatically by the asymptotic normal value. For a two-sided central 95\% interval, the total tail probability is $0.05$, which is split equally between both tails of the symmetric Student-$t$ distribution. Therefore, the corresponding coverage factor is the $97.5$th percentile, $t_{0.975,\nu}$. Thus, for a Type A standard uncertainty associated with a mean quantity, the 95\% expanded uncertainty is written as:

$$
U_{95}=t_{0.975,\nu}\,u_A(\overline{x}).
$$

As $\nu\rightarrow\infty$, this factor tends to the normal-distribution value $1.96$, which explains why the approximation $k\approx2$ becomes acceptable for large samples.

In the present context, however, the basic observations are not usually independent repetitions of a scalar experiment, but finite time series that may exhibit temporal correlation. For a mean quantity estimated from such a record, the raw sample size $N$ overestimates the actual amount of independent information. A practical and physically meaningful approximation is therefore to define the degrees of freedom of the Type A contribution from the effective sample size:

$$
\nu_{x,\mathrm{eff}}\approx N_{\mathrm{eff},x}-1.
$$

Accordingly, the 95\% random interval associated with the mean of a correlated time series may be written as:

$$
U_{95,A}(\overline{x})=t_{0.975,\nu_{x,\mathrm{eff}}}\frac{s_x}{\sqrt{N_{\mathrm{eff},x}}}.
$$

This expression is the natural extension of the usual standard error of the mean to the case of mildly correlated laboratory signals.

(welch_satterthwaite_formula)=
### The Welch$-$Satterthwaite formula

For a generic measurand $y=f(x_1,\dots,x_n)$, the final combined standard uncertainty is typically obtained from several input contributions, each one weighted by the corresponding sensitivity coefficient:

$$
u_c^2(y)=\sum_{i=1}^{n} c_i^2u_i^2,
\qquad
c_i=\frac{\partial f}{\partial x_i}.
$$

If those input uncertainty components do not all possess the same degrees of freedom, the effective number of degrees of freedom associated with the combined standard uncertainty is approximated through the **Welch$-$Satterthwaite formula**:

$$
\nu_{\mathrm{eff}}=
\frac{u_c^4(y)}
{\displaystyle\sum_{i=1}^{n}\frac{c_i^4u_i^4}{\nu_i}}.
$$

This effective value is then used to select the coverage factor for the desired reporting level. In particular, the expanded uncertainty associated with a 95\% interval is written as:

$$
U_{95}(y)=t_{0.975,\nu_{\mathrm{eff}}}\,u_c(y).
$$

This is the rigorous route by which the widely used 95\% interval is incorporated into the linear propagation framework. In large-sample situations, the factor approaches 2; in small-sample or low-degree-of-freedom cases, it may be significantly larger, and replacing it blindly by 2 may lead to underestimation of the reported interval.

For many Type B components, the associated degrees of freedom may be treated as very large, or effectively infinite, when the corresponding standard uncertainties are regarded as known with negligible uncertainty relative to the Type A terms. Under those conditions, the Welch$-$Satterthwaite formula is mainly controlled by the finite-degree-of-freedom components, which in the present framework are typically the Type A contributions extracted from time records.

(note_on_welch_satterthwaite)=
:::{note}
The **Welch$-$Satterthwaite formula** is used to assign an effective number of degrees of freedom to a **combined standard uncertainty** that is itself obtained by combining several uncertainty contributions with different degrees of freedom.

To see why the formula has the form it does, consider the combined standard uncertainty of a measurand $y$:

$$
u_c^2(y)=\sum_{i=1}^{n} c_i^2u_i^2,
$$

where:

$$
c_i=\frac{\partial f}{\partial x_i}
$$

is the sensitivity coefficient associated with the $i$-th input quantity, and $u_i$ is the corresponding standard uncertainty contribution.

Define:

$$
q_i=c_i^2u_i^2.
$$

Then the propagated variance estimate can be written as:

$$
u_c^2(y)=\sum_{i=1}^{n} q_i.
$$

The difficulty is that each $q_i$ may itself be an estimated variance contribution, and therefore may be associated with a different number of degrees of freedom $\nu_i$. The sum of such terms does not, in general, follow an exact scaled chi-squared distribution. The Welch$-$Satterthwaite idea is to replace that sum by a **single equivalent variance estimate** having the same mean and approximately the same variance, and then assign to that equivalent quantity an effective number of degrees of freedom $\nu_{\mathrm{eff}}$.

Let the surrogate variance estimate be denoted by $Q^\ast$, and suppose it behaves like a scaled chi-squared variable with effective degrees of freedom $\nu_{\mathrm{eff}}$. For such a quantity, the relative variance is approximately:

$$
\frac{\operatorname{Var}(Q^\ast)}{\left[\mathbb E(Q^\ast)\right]^2}\approx\frac{2}{\nu_{\mathrm{eff}}}.
$$

Likewise, for each individual component $q_i$ associated with $\nu_i$ degrees of freedom, one has the analogous approximation:

$$
\frac{\operatorname{Var}(q_i)}{\left[\mathbb E(q_i)\right]^2}\approx\frac{2}{\nu_i}.
$$

Assuming the contributions are independent, the variance of the sum is:

$$
\operatorname{Var}\left(u_c^2(y)\right)=\sum_{i=1}^{n}\operatorname{Var}(q_i)\approx2\sum_{i=1}^{n}\frac{q_i^2}{\nu_i}.
$$

Now impose the defining condition of the equivalent approximation: the surrogate $Q^\ast$ must have the same mean as the actual combined variance estimate,

$$
\mathbb E(Q^\ast)=u_c^2(y),
$$

and the same approximate variance,

$$
\operatorname{Var}(Q^\ast)\approx\operatorname{Var}\!\left(u_c^2(y)\right).
$$

Using the scaled chi-squared form for $Q^\ast$ gives:

$$
\operatorname{Var}(Q^\ast)\approx\frac{2u_c^4(y)}{\nu_{\mathrm{eff}}}.
$$

Equating both variance expressions:

$$
\frac{2u_c^4(y)}{\nu_{\mathrm{eff}}}=2\sum_{i=1}^{n}\frac{q_i^2}{\nu_i}.
$$

After cancelling the factor 2 and substituting back $q_i=c_i^2u_i^2$, one obtains:

$$
\nu_{\mathrm{eff}}=
\frac{u_c^4(y)}
{\displaystyle\sum_{i=1}^{n}\frac{c_i^4u_i^4}{\nu_i}}.
$$

Thus, the structure of the formula reflects a moment-matching argument:
- the numerator involves $u_c^4(y)$ because the comparison is performed at the level of the **variance of a variance estimate**,
- the denominator weights each contribution according to both its propagated magnitude ($c_i^2u_i^2$) and the uncertainty in that estimate (through $\nu_i$),
- and components with small $\nu_i$ penalize $\nu_{\mathrm{eff}}$ more strongly, because their uncertainty estimates are themselves less reliable.

Once $\nu_{\mathrm{eff}}$ has been obtained, it is used to determine the appropriate Student-$t$ coverage factor for the desired reporting level. In particular, the expanded uncertainty associated with a two-sided 95\% interval is written as:

$$
U_{95}(y)=t_{0.975,\nu_{\mathrm{eff}}}u_c(y).
$$

Thus, the Welch$-$Satterthwaite formula provides the bridge between the combined standard uncertainty and the final 95\% reported interval whenever the constituent uncertainty contributions do not all possess the same degrees of freedom.
:::

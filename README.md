# Two-layer TGS fitting

This MATLAB workflow fits the film diffusivity `alpha_f`, substrate
diffusivity `alpha_s`, and interface resistance `R` for one measured spot at
a time.

All three physical parameters are mandatory in every fit. They cannot be
enabled or disabled individually; `S.p0`, `S.pLower`, and `S.pUpper` must
each contain exactly `[alpha_f, alpha_s, R]`.

Run `MainTwoLayerScript` after editing its user-input block. `S.calFile` must
name exactly one calibration file, and the `run_name` entries in that file
select the one spot to fit. Other spots in `S.dataDir` are ignored. To fit a
different spot, manually change `S.dataDir` as needed and select that spot's
calibration file.

The data loader pairs POS and NEG files by nominal grating, spot, and run.
Every complete run is retained as an independent observation; repeated runs
are not averaged before fitting.

One complete run is sufficient to fit all three material parameters because
it contains many intensity-versus-time samples. Each additional run is fit
independently and then used to refine the spot-level result.

The fitting sequence follows the single-layer TGS workflow:

- fit a smoothed thermal response to initialize the material parameters;
- fit the unsmoothed intensity trace with the bounded two-layer temperature,
  displacement, and offset model;
- repeat that thermal-model fit independently for each run.

Temperature, displacement, and offset amplitudes are profiled independently
inside each run's fit, leaving only `alpha_f`, `alpha_s`, and `R` in the
nonlinear optimizer. The workflow accepts one or more runs from exactly one
spot and one grating.

The model contains no SAW frequency, damping time, or oscillatory fitting
term. Any acoustic oscillation present in a measured trace remains in the
reported residual and can increase the fitted-parameter uncertainty.

The reported spot-level value is the arithmetic mean of the identifiable
run estimates. Following the single-layer `meanAndError` convention, the
reported uncertainty combines the propagated within-fit errors with the
between-run standard deviation. Results include a per-run table so every
individual fit remains inspectable. A run remains in that table even when a
bound or rank diagnostic excludes one of its parameter estimates from the
spot-level mean.

The output reports bound activity, Jacobian rank and condition, local
sensitivity, and standard errors. A parameter at a bound is marked not
identifiable and its symmetric standard error is reported as `NaN`.

Run the regression tests with:

```matlab
results = runtests("test_single_grating_repeated_runs.m");
assertSuccess(results)
```

# LCA equations

The functional unit is treatment of 1 t dry mixed feedstock consisting of
0.5 t dry sewage sludge and 0.5 t dry rice straw.

For hydrothermal pathways, net GWP is calculated as:

```text
net GWP = transport and pretreatment
        + electricity
        + process-gas emissions
        + methane-leakage emissions
        + supplementary natural-gas heat emissions
        - avoided mineral fertilizer production
        - soil carbon sequestration
```

The soil-carbon credit is:

```text
soil-carbon credit = - hydrochar yield
                     x measured C fraction
                     x 44/12
                     x carbon-stability factor
```

For direct incineration:

```text
net GWP = transport
        + direct incineration emissions
        - exported electricity credit
```

One-at-a-time sensitivity analysis fixes all parameters at EFS values and
varies one parameter at a time between its lower- and higher-emission endpoint.
The sensitivity range is the absolute difference between the two endpoint net
GWP values.

Break-even heat recovery is the minimum heat-recovery efficiency at which PHC
net GWP equals zero for a specified PHC carbon-stability factor, with other
parameters fixed at EFS values.

Monte Carlo analysis uses triangular uncertainty ranges bounded by CSS and OIS
values, with EFS as the mode. The reported 95% interval crosses zero
(-264.0 to +24.0 kg CO2-eq t-1), so net-negative operation is likely but not
guaranteed.


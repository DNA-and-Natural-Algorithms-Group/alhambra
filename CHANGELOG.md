- MultiFixedConcentration allows a `min_volume` setting, which will raise an error of the minimum
  volume to be transferred for any component is too low.
- Volume and concentration settings can now be changed, not just initialized, as strings.
- Reference data for mixes is now its own class, `Reference`, and rounds to 1e-6 nM.
- Mixes now use Decimal instead of floats throughout for units, solving floating point errors.

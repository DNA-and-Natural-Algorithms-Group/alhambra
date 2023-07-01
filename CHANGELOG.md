## a16

- Add a seed parameter to set a single seed for a tileset.

## a15

- Allow a default seed offset for simulations to be included in seed specifications (if
  implemented by the seed class).
- Preliminary support for DX tiles and older files.

## a14

- Internal SVG drawing routines, removing drawSvg_svgy and lxml dependencies.
- Xgrow release dependency.
- Support for stickydesign-calculated sequence-dependent glues for xgrow/rgrow export (fixes regression from 1.0).
- More SST motifs.
- More flatish motifs (seeds, extensions on single tiles, etc).
- Early support for rgrow export.
- Improved error messages.

## a13

- Python 3.11 compatibility.

## a12

- Mix-to-tileset code moved to TileSet.

## Previous

- MultiFixedConcentration allows a `min_volume` setting, which will raise an error of the minimum
  volume to be transferred for any component is too low.
- Volume and concentration settings can now be changed, not just initialized, as strings.
- Reference data for mixes is now its own class, `Reference`, and rounds to 1e-6 nM.
- Mixes now use Decimal instead of floats throughout for units, solving floating point errors.
- Formatting for mix table entries now takes place in `MixLine`.


# v1.1.0

Fixes broken workaround for ruamel.yaml bug now fixed upstream, adds double-tile sensitivity, includes
seed file to make xor example work.

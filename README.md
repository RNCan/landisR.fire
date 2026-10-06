# landisR.fire

Development package for modular fire models that interact with landisR.

The first scientific target is strict reproduction of LANDIS-II Base Fire.
Base Fire will be implemented as one configuration of replaceable modules for ignition, target size, spread, severity and mortality.

The architecture is intentionally open so users can later substitute empirical fire-size distributions, composition-sensitive burning/severity, alternative spread algorithms or completely custom functions without modifying landisR.

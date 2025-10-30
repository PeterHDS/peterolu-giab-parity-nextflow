# Portable Precision — GIAB Parity

Scaffold created.

### Lessons learned (Peter Olu DS)
- Channels in DSL2 are single-consumer. If two steps need the same data, either split properly or give each branch its own input.
- If a module is used twice, **alias it** (e.g. `GATK_HC_A` and `GATK_HC_B`).
- Make reference/index steps **idempotent**; don’t publish intermediates that can collide.
- hap.py needs a real output folder and clean, indexed VCFs with matching contigs and headers.
- Pin container versions and stick to LF line endings for predictable runs.

For the full troubleshooting story, see: `docs/troubleshooting.md`.

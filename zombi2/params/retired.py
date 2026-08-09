"""The spellings that were the grammar and are not any more, with the sentence each one needs.

A user porting a script hits the old vocabulary in two places: in **text** (a ``--birth`` flag, a
``--params`` file) and in **Python** (an import, a keyword argument). Both have to answer with the
replacement rather than with "unknown name" or "unexpected keyword argument", which name the mistake
without naming the fix — so the table lives here, on its own, and `parse` and the package's
``__getattr__`` both read it. One table, two surfaces: a retirement recorded once cannot be answered
well in one place and badly in the other.

A difflib guess is no help for any of these: the replacement is usually a **verb** on the parameter
rather than another name, so there is nothing close enough to guess at.
"""

from __future__ import annotations

#: Names that were the written form and are not any more, each with what to write instead.
RETIRED = {
    "DrivenBy": ("write the verb that says what the number does: .scaled_by(driver, mapping) on a "
                 "rate or an extent, .weighted_by(driver, mapping) on a Recipients() rule, "
                 ".set_by(driver, mapping) to replace the base rather than scale it"),
    "ScaledBy": ("the verbs are methods on the parameter now: "
                 "PerCopy(0.25).scaled_by(driver, mapping)"),
    "SetBy": ("the verbs are methods on the parameter now, and a replaced base is written from the "
              "bare scope: PerCopy().set_by(driver, mapping)"),
    "Weights": ("a choice is written from its own entry point: "
                "Recipients().weighted_by(driver, mapping)"),
    "OnTime": ("the run's clock has its own verb: PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})"),
    "OnTotalDiversity": ("write the driver and scale by it: "
                         "PerLineage(1.0).scaled_by(TotalDiversity(cap=100))"),
    "Drawn": ("a random value has its own verb, and its unit goes plural: "
              "PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))"),
    "Inherited": ("a random value has its own verb, and inheritance is a law of it: "
                  "PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))"),
    "ByFamily": ("write PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))"),
    "ByLineage": ("write PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))"),
    "FromParent": ("write PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))"),
}

#: Keywords that were the written form and are not any more. They are caught here rather than by the
#: class or the verb they were passed to, because "unexpected keyword argument" names neither what
#: the argument meant nor what replaced it.
RETIRED_KEYWORDS = {
    "spread": ("write the law out, so the page says which distribution it is and what it "
               "distributes — varying_among('families', LogNormal(0.0, spread)) for the value, "
               "varying_among('lineages', Drift(LogNormal(0.0, spread))) for the per-split step. "
               "One word named both, which is why it is gone"),
    "per": ("the unit is varying_among's first argument, and units are plural: "
            "varying_among('families', LogNormal(0.0, 0.5))"),
}


def name_message(name: str) -> str:
    """The sentence for a retired **name**, wherever it is written."""
    return f"{name} is no longer a rate name — {RETIRED[name]}"


def keyword_message(keyword: str) -> str:
    """The sentence for a retired **keyword argument**, wherever it is passed."""
    return f"{keyword}= is no longer written — {RETIRED_KEYWORDS[keyword]}"


def check_no_retired_keywords(kwargs: dict, *, where: str) -> None:
    """Refuse a retired keyword by name, and refuse every other leftover keyword as Python would.

    A verb that takes ``**retired`` so it can answer ``per=`` and ``spread=`` properly would
    otherwise **swallow a typo**: ``varying_among('families', LogNormal(0.0, 0.5), bns=8)`` would
    run with no binning and no complaint. So nothing is allowed through — a retired keyword gets the
    sentence naming its replacement, anything else gets the message Python would have given.
    """
    for name in kwargs:
        if name in RETIRED_KEYWORDS:
            raise TypeError(keyword_message(name))
    if kwargs:
        raise TypeError(
            f"{where}() got an unexpected keyword argument {sorted(kwargs)[0]!r}")


__all__ = ["RETIRED", "RETIRED_KEYWORDS", "name_message", "keyword_message",
           "check_no_retired_keywords"]

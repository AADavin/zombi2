# Matching empirical profiles (ABC)

!!! warning "Experimental — withheld from this release"
    ABC profile-matching inference (`match_profiles`, `match_profiles_smc`, `match_coupled`
    and the summary statistics) is **not part of the public v1 API** yet: it is not exported
    from the top-level `zombi2` namespace, has no CLI command, and its reference docs are
    withheld pending stabilisation and documentation.

    The implementation still ships in-tree and is fully tested. If you want to experiment with
    it, import it explicitly:

    ```python
    from zombi2.matching import match_profiles, match_profiles_smc, cooccurrence_summary
    ```

    The API may change without notice. This page will be restored when the module is promoted
    to the stable surface.

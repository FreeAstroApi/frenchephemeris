# Contributing to FrenchEphemeris

Thank you for helping. Please read this page before opening a pull request.

## Why there is a contributor licence agreement

FrenchEphemeris is **dual-licensed**: AGPL-3.0-only for everyone, and a paid commercial
licence sold by the copyright holder, Gabriel Naulleau EI, trading as FreeAstroAPI (see
[COMMERCIAL-LICENSE.md](COMMERCIAL-LICENSE.md)). The copyright holder can only offer the
commercial licence for code it has the right to license. A contribution therefore needs
the contributor's agreement below before it can be merged. Pull requests without it
cannot be accepted.

## How to agree

Add this line to the description of your first pull request, and keep the `Signed-off-by`
trailer on your commits (`git commit -s`):

```
I have read the FrenchEphemeris Contributor Licence Agreement in CONTRIBUTING.md
(version 1.0) and I agree to it for this and all my future contributions.
```

When the CLA assistant is enabled on the repository, it asks you to sign once, records
your signature and checks later pull requests automatically.

Contributions made on behalf of an employer, or as part of your job, need your employer's
agreement as well: write to contact@freeastroapi.com first.

Please also follow the [code of conduct](CODE_OF_CONDUCT.md), and report security issues
privately as described in [SECURITY.md](SECURITY.md).

## Technical rules

- Build: `cargo build --workspace`; tests: `cargo test --release --workspace --exclude fe-py`
  and `pytest tests` (packs: `python -m frenchephemeris.packs download`, see README.md). The public test suite must pass **without** the private
  Swiss-oracle corpus (`private/oracle/`, `FE_ORACLE_DIR`).
- **Never** add Swiss Ephemeris source code, data files or outputs to the public tree.
  Fixtures that contain Swiss outputs are written to the private corpus (see
  `tools/privatize_oracle_fixtures.py`; run it with `--check` before you push).
- Every numeric constant must come from a published source or from our own fits to
  independent data (JPL, ERFA, IERS, NASA, IAU); record the source in the module docs and
  the relevant `docs/*REPORT*.md` provenance table.
- New source files start with the SPDX header:
  ```
  // SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-FrenchEphemeris-Commercial
  // Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI
  ```
  (`#` instead of `//` in Python). Generated tables add a `// Data:` line naming their source.
- Third-party code or data: only with a licence compatible with **both** the AGPL-3.0 and
  proprietary redistribution (e.g. MIT, BSD, Apache-2.0, public domain, CC BY), and listed
  in [NOTICE](NOTICE).

---


## FrenchEphemeris Contributor Licence Agreement — version 1.0 (30 October 2026)

This agreement is between **you** (the person, or the entity on whose behalf you act, that
submits a Contribution) and **Gabriel Naulleau, entrepreneur individuel, trading as
FreeAstroAPI**, SIRET 917 923 096 00015, 48 rue du Montparnasse, 75014 Paris, France, and
any successor to whom the project is transferred under section 8 (**"we"**).

### 1. Definitions

- **"Project"**: FrenchEphemeris, as published by us at
  https://github.com/FreeAstroApi/frenchephemeris or any successor location.
- **"Contribution"**: any code, data, documentation or other material that you submit to us
  for inclusion in the Project, by pull request, patch, issue or any other means, with this
  agreement in force. Material you mark "Not a Contribution" when you submit it is excluded.

### 2. Copyright licence

2.1 You grant us, for each Contribution, a non-exclusive, royalty-free, worldwide licence,
for the full legal term of copyright protection (including any extension), of the
following rights:

  a. the right to reproduce the Contribution, permanently or temporarily, in whole or in
     part, by any means and in any form;
  b. the right to translate, adapt, arrange and modify it, and to reproduce the resulting
     works;
  c. the right to distribute it and to communicate it to the public, including by making it
     available online, by sale, rental, free distribution or as a hosted service;
  d. the right to include it, modified or not, in other works;

for the purpose of developing, publishing, maintaining and exploiting the Project and works
based on it, under any licence terms we choose, including the GNU AGPL-3.0, other
open-source licences and proprietary commercial licences.

2.2 We may sublicense these rights, through several tiers, to licensees of the Project, and
transfer them under section 8.

2.3 You grant this licence free of charge and acknowledge that no remuneration is due for
it.

2.4 You keep the copyright in your Contributions and may use them for any purpose,
including under other licences.

### 3. Patent licence

You grant us and every recipient of the Project a non-exclusive, royalty-free, worldwide
licence, for the life of the patents concerned, under any patent claims you can license
that are necessarily infringed by your Contribution alone or by its combination with the
Project as it stood when you submitted it, to make, use, sell, offer for sale, import and
otherwise transfer the Project. If a recipient starts litigation claiming that the Project
or a Contribution infringes a patent, the licences granted to that recipient under this
section end on the date the litigation is filed.

### 4. Our commitment

Every Contribution that we include in a public release of the Project remains available to
the public under the GNU AGPL-3.0 or another licence approved by the Open Source Initiative,
in addition to any commercial licence we offer.

### 5. Moral rights

French law gives authors moral rights that cannot be waived. For software, those rights are
limited by article L121-7 of the Code de la propriété intellectuelle: you may not object to
modifications or adaptations of your Contribution made under section 2 unless they harm
your honour or reputation, and you may not withdraw it. You agree that you may be credited
in the Project's history and contributor list in the usual way of open-source projects, and
that we are not required to name you elsewhere.

### 6. Your statements

You confirm that:

  a. each Contribution is your original creation, or you have the right to submit it under
     this agreement; any third-party material it contains is identified in your submission,
     with its source and licence, and is available under a licence compatible with both the
     AGPL-3.0 and proprietary redistribution;
  b. it contains no code, data or output of the Swiss Ephemeris, and was not written by
     reading Swiss Ephemeris source code;
  c. if you are employed or under contract, your employer or client has authorised your
     Contribution or has no rights in it;
  d. you are at least 18 years old, or the holder of parental authority has agreed to this
     agreement for you;
  e. you will tell us promptly if any of these statements becomes inaccurate.

Except for these statements, Contributions are provided "as is", without warranty.

### 7. No obligation

We are not required to include a Contribution in the Project, or to keep it there.

### 8. Transfer

We may transfer the rights granted to us under this agreement to a company we control or
that succeeds to the Project (for example on incorporation of our sole proprietorship),
which then takes over our commitments in section 4.

### 9. Signature and records

You accept this agreement by the statement in your pull request, by signing with the CLA
assistant, or by a signed copy sent to contact@freeastroapi.com; each is an electronic
signature you agree to be bound by. We keep a record of your name, GitHub account, email
address, the version accepted and the date, as controller under the GDPR (Regulation (EU)
2016/679), to prove the rights granted, for as long as the Project contains your
Contributions and the legal limitation periods after that. You may exercise your data
rights by writing to contact@freeastroapi.com and may complain to the CNIL.

### 10. Versions

This agreement applies to all your Contributions submitted while it is in force. A later
version applies only to Contributions submitted after you accept it.

### 11. Governing law

This agreement is governed by French law. The courts within the jurisdiction of the Paris
Court of Appeal (cour d'appel de Paris) have exclusive jurisdiction, except where
mandatory rules give you the right to sue or be sued in the courts of your residence.

# FrenchEphemeris commercial licence

SPDX identifier of this licence: `LicenseRef-FrenchEphemeris-Commercial`

This page has two parts: a plain-language overview, and the **FrenchEphemeris Commercial
Licence Agreement, version 1.0**, which is the binding text. If the overview and the
Agreement differ, the Agreement prevails.

## Licensor

Gabriel Naulleau, entrepreneur individuel (EI, French sole proprietorship,
micro-entreprise), trading as **FreeAstroAPI**. "FreeAstroAPI" is a trade name, not a
separate company.

- SIRET: 917 923 096 00015
- Registered address: 48 rue du Montparnasse, 75014 Paris, France
- VAT: TVA non applicable, art. 293 B du CGI (VAT franchise)
- Contact: [contact@freeastroapi.com](mailto:contact@freeastroapi.com)

Copyright (C) 2026 Gabriel Naulleau EI, trading as FreeAstroAPI.

## Overview

FrenchEphemeris is available under **either**:

1. the **GNU Affero General Public License, version 3 only** (`AGPL-3.0-only`, see
   [LICENSE](LICENSE)), free of charge; **or**
2. this **commercial licence**, for USD 200, one-time, per legal entity.

The commercial licence lets you use FrenchEphemeris without the AGPL's copyleft
obligations. It does not change the third-party terms listed in [NOTICE](NOTICE) (JPL,
ERFA BSD-3 notice, IAU WGSN CC BY 4.0 attribution, ESA/CDS acknowledgements, Rust crate
licences), which apply under both options.

### What the AGPL requires (summary, not a substitute for the licence text)

- **Distribution.** If you convey FrenchEphemeris or a work based on it (an app, a library,
  a container image, a Python wheel), you must license the whole work under the AGPL-3.0
  and provide its complete corresponding source code (AGPL §§ 4–6).
- **Network use (§ 13).** If you modify FrenchEphemeris and let users interact with it
  remotely (a web API, a SaaS, a website computing charts on your server), you must offer
  those users the complete corresponding source code of your modified version.
- **Notices.** Keep the copyright, licence and NOTICE files; mark your changes.

Research, teaching and open-source use are free under the AGPL-3.0 when its conditions
are met. The AGPL does not forbid commercial use; the commercial licence is for those who
do not want to publish their own source code under the AGPL.

### Who needs the commercial licence

You probably need it if you use FrenchEphemeris **and** do not want to meet the AGPL
obligations, for example in:

- a closed-source application, plug-in or service that ships FrenchEphemeris (desktop,
  mobile, embedded);
- a hosted service or API built on a modified FrenchEphemeris whose source you do not want
  to publish under the AGPL.

You do **not** need it if you comply with the AGPL-3.0. If you are unsure, write to
contact@freeastroapi.com before you ship.

### Business terms at a glance

| | Commercial licence |
|---|---|
| Customers | Businesses and professionals only (B2B); not offered to consumers |
| Price | USD 200, one-time, excluding taxes |
| Duration | The full legal term of copyright in the Software |
| Updates | All versions the licensor releases under the FrenchEphemeris name, at no extra cost; no support or service-level obligation |
| Scope | One legal entity: unlimited products, servers, deployments and developers, including contractors working solely on the licensee's products |
| Affiliates | Not included; each company needs its own licence (group licences on request) |
| Downstream use | Unlimited customers and end users of your products, who need no licence of their own |
| Not included | OEM use (see the test below): needs a separate OEM licence, priced on request |
| Proof of licence | A licence certificate naming the licensee, date and licence number |

### OEM test

The standard licence covers *products built with* FrenchEphemeris. It does **not** cover
offerings whose primary value to the recipient is FrenchEphemeris itself. An offering is
**OEM use** if any of these is true:

1. it is a **library, SDK, plug-in, component or source distribution** that other
   developers integrate into their own software;
2. it is an **API or data feed whose primary output is raw astronomical data** (planetary
   or point positions, speeds, house cusps, ayanamsas, event times) offered to third
   parties to build their own astrology or astronomy products;
3. it is a **white-label or rebranded** version of FrenchEphemeris or of its outputs.

Covered by the standard licence: a horoscope or chart-reading API whose responses include
positions as part of interpretations or charts; a website, SaaS or app for end users;
internal tools. Needs an OEM licence: an "ephemeris API" selling raw positions per call to
other developers; a package wrapping FrenchEphemeris published for others; an astrology
engine sold to other software vendors.

### How to buy

Email [contact@freeastroapi.com](mailto:contact@freeastroapi.com) with your legal entity
name, registered address, country, VAT or company identification number, and one line
describing the product that will include FrenchEphemeris. You receive an order form
referring to this Agreement and an invoice; the licence starts when the order form is
accepted and the price paid, and you then receive the licence certificate.

---

# FrenchEphemeris Commercial Licence Agreement — version 1.0 (30 October 2026)

**Between** Gabriel Naulleau, entrepreneur individuel, trading as FreeAstroAPI, SIRET
917 923 096 00015, 48 rue du Montparnasse, 75014 Paris, France (the **"Licensor"**), **and**
the legal entity named in the Order Form (the **"Licensee"**).

## 1. Definitions

- **"Agreement"**: this text, the Order Form and the Licence Certificate.
- **"Order Form"**: the document (quote, order form or written order confirmation) that
  names the Licensee, refers to this Agreement and states the price.
- **"Licence Certificate"**: the document issued by the Licensor after payment, stating
  the Licensee, the licence number and the Effective Date.
- **"Effective Date"**: the date on which the Licensor receives full payment of the price.
- **"Software"**: FrenchEphemeris in source and object code, its documentation and the
  coefficient packs distributed by the Licensor, in every version released by the Licensor
  under the FrenchEphemeris name before or after the Effective Date.
- **"Licensed Product"**: a product or service of the Licensee that incorporates or uses
  the Software, other than OEM Use, offered to the Licensee's customers or end users or
  used internally.
- **"End User"**: any person or entity that uses a Licensed Product.
- **"OEM Use"**: an offering that meets the OEM test in the overview above (libraries, SDKs
  and source distributions for other developers; APIs or data feeds whose primary output is
  raw astronomical data; white-label versions).
- **"Contractor"**: a person or company working for the Licensee solely on the Licensee's
  Licensed Products.

## 2. Formation and acceptance

2.1 The Agreement is formed when the Licensee accepts the Order Form, by signature
(handwritten or electronic) or by written confirmation by email, and takes effect on the
Effective Date. Before accepting, the Licensee has had this Agreement available to read and
keep at the address given in the Order Form.

2.2 The Licensee states that it acts for the purposes of its trade, business or profession
and that the person accepting has authority to bind it. The Software is not licensed to
consumers under this Agreement.

2.3 This Agreement replaces any purchase conditions of the Licensee, even if they are
referred to in a purchase order, unless the Licensor accepts them in writing.

## 3. Grant of rights

3.1 From the Effective Date, the Licensor grants the Licensee a worldwide, non-exclusive
licence, for the full legal term of copyright protection in the Software (including any
extension), to:

  a. reproduce the Software, permanently or temporarily, in whole or in part, by any means
     and in any form, for loading, displaying, running, transmitting and storing it;
  b. translate, adapt, arrange and modify the Software and reproduce the resulting
     software;
  c. integrate the Software, modified or not, into Licensed Products;
  d. distribute the Software in object code form, and in source form where needed for
     interpreted languages, only as an integrated part of Licensed Products, to customers
     and End Users, in any number, by sale, subscription, free distribution or as a hosted
     service;
  e. let Contractors exercise these rights for the Licensee.

3.2 These rights are granted for the purpose of developing, operating, marketing and
supporting Licensed Products, outside the terms of the AGPL-3.0.

3.3 End Users may use the Software as part of Licensed Products without a licence of their
own. The Licensee may grant them the rights needed for that use, under terms at least as
protective of the Licensor as sections 5 and 8, but may not grant them a licence to the
Software on its own.

3.4 The Licensor keeps all rights not expressly granted. Nothing in this Agreement transfers
ownership of the Software or grants rights to the Licensor's trade names, which are covered
by [TRADEMARKS.md](TRADEMARKS.md).

## 4. Updates and support

4.1 The licence covers every version of the Software the Licensor releases under the
FrenchEphemeris name, before or after the Effective Date, at no extra cost.

4.2 The Licensor has no obligation to release updates, to maintain any feature, or to
provide support, training or any service level. Security vulnerabilities are handled as
described in [SECURITY.md](SECURITY.md).

## 5. Restrictions

The Licensee shall not:

  a. engage in OEM Use without a separate OEM licence;
  b. remove or alter the copyright, licence, NOTICE or attribution statements in the
     Software;
  c. present the Software, or a product based on it, as made or endorsed by the Licensor,
     except as allowed by TRADEMARKS.md;
  d. sublicense, rent or transfer the Software except as allowed by sections 3 and 13.

The Software may also be used under the AGPL-3.0; nothing in this Agreement restricts rights
the Licensee holds under that licence.

## 6. Price, invoicing and payment

6.1 The price is stated in the Order Form (USD 200 at the date of this version). It is a
one-time fee, payable in full before the Effective Date, by bank transfer or the payment
method stated on the invoice. It is non-refundable except under sections 8.2 and 10.4.

6.2 Prices exclude taxes. The Licensor applies the French VAT franchise (art. 293 B of the
Code général des impôts) and charges no French VAT. Where the Licensee is a business
established outside France, any VAT or similar tax due in the Licensee's country is
accounted for by the Licensee under its local rules (in the European Union, under the
reverse-charge mechanism). Any withholding tax is borne by the Licensee, so that the
Licensor receives the full price.

6.3 Invoices are issued in accordance with French law, electronically where required. If
an amount is paid late, the Licensee owes, without reminder, late-payment interest at three
times the French statutory interest rate (art. L441-10 of the Code de commerce) and the
fixed recovery indemnity of EUR 40 (art. D441-5 of the Code de commerce), without prejudice
to further compensation for documented recovery costs. Until full payment, no licence is
granted.

## 7. Third-party components and data

7.1 The Software contains or uses third-party components and data listed in
[NOTICE](NOTICE), under their own terms, which continue to apply.

7.2 With every distribution of a Licensed Product containing the Software, the Licensee
reproduces the NOTICE attributions (including the ERFA BSD-3 notice and the IAU WGSN
CC BY 4.0 credit), for example in an "About" or "Licences" screen or file.

## 8. Warranties

8.1 The Licensor warrants that it has the right to grant the licence in section 3.

8.2 If a third party claims that the Software, as delivered by the Licensor, infringes its
intellectual property rights, the Licensee notifies the Licensor promptly and the Licensor
may, at its choice and expense: obtain the right for the Licensee to continue using the
Software; modify or replace the affected part so that it no longer infringes, without
material loss of function; or, if neither is reasonably possible, terminate the licence for
the affected part and refund the price paid. This is the Licensee's sole remedy for such
claims, subject to section 9.2. It does not apply to claims arising from modifications made
by or for the Licensee, from combinations with other software or data, or from use in
breach of this Agreement.

8.3 Apart from section 8.1, the Software is provided "as is". To the extent permitted by
law between professionals, the Licensor gives no other express or implied warranty,
including of accuracy, uninterrupted or error-free operation, or fitness for a particular
purpose. The Licensee is responsible for choosing the Software, for testing it in its own
products, and for the results it obtains. The Software is not designed for navigation,
surveying, safety-critical, medical or legal timekeeping uses.

## 9. Liability

9.1 Each party is liable only for direct and foreseeable damage caused by its fault. The
Licensor is not liable for loss of profit, revenue, business, data or goodwill, or for
damage caused by the Licensee's products, by third-party components, or by use of the
Software outside its documentation.

9.2 The Licensor's total liability under or in connection with this Agreement is limited
to the price paid by the Licensee. This limitation does not apply to liability for fraud,
wilful misconduct or gross negligence, for death or personal injury, under the rules on
liability for defective products, or to any other liability that cannot be limited or
excluded under applicable law. The parties agree that this limitation reflects the
one-time price and the allocation of risk in this Agreement and does not deprive the
Licensor's essential obligation of its substance.

9.3 Any claim must be made within one year from the date the claiming party became aware,
or should have become aware, of the facts giving rise to it, to the extent permitted by
law.

## 10. Term and termination

10.1 The Agreement runs from the Effective Date for the term stated in section 3.1.

10.2 Either party may terminate the Agreement by written notice if the other party commits
a material breach and does not remedy it within 30 days of a written notice describing the
breach. Use of the Software outside the scope of section 3 or in breach of section 5 is a
material breach.

10.3 On termination for the Licensee's breach, the Licensee stops distributing the Software
and Licensed Products containing it and stops using it in hosted services within 30 days.
Copies already distributed to End Users before termination may continue to be used by
them.

10.4 If the Licensee terminates for the Licensor's uncured material breach, the Licensor
refunds the price paid.

10.5 Sections 5, 7, 8.3, 9, 10.3 to 10.5, 14 and 15 survive termination.

## 11. Security and regulatory compliance

11.1 The Licensor handles vulnerabilities in the Software as described in SECURITY.md,
including coordinated disclosure and, where required, notification of authorities under
the EU Cyber Resilience Act (Regulation (EU) 2024/2847). Security fixes are released as
updates covered by section 4.

11.2 The Licensee remains solely responsible for the compliance of its own Licensed
Products, including its obligations as their manufacturer, provider or operator, and for
applying updates.

11.3 Each party complies with the export-control and sanctions laws that apply to it. The
Licensee states that it, its owners and its End Users it knowingly serves are not subject to
European Union, United Nations, French, United Kingdom or United States sanctions that would
prohibit this Agreement. The Licensor may refuse an order, or suspend the licence, to the
extent needed to comply with such laws.

## 12. Personal data

12.1 The Licensor processes the contact and billing details of the Licensee's
representatives as controller, to perform this Agreement, issue invoices and meet its legal
(accounting and tax) obligations, under the GDPR (Regulation (EU) 2016/679). The data are
kept for the duration of the Agreement and then for the periods required by law (ten years
for accounting records). Data subjects may exercise their rights by writing to
contact@freeastroapi.com and may lodge a complaint with the CNIL.

12.2 The Software runs on the Licensee's systems; the Licensor does not receive or process
data processed by the Licensee with the Software, and is not its processor.

## 13. Transfer

13.1 The Licensee may not transfer this Agreement, except to a successor that acquires all
or substantially all of its business, or the Licensed Product that uses the Software, and
that accepts this Agreement in writing; the Licensee notifies the Licensor of the transfer.

13.2 The Licensor may transfer this Agreement, together with the rights in the Software, to
a company it controls or that succeeds to its business (for example on incorporation of its
sole proprietorship). The Licensee's rights continue unchanged.

## 14. General

14.1 **Entire agreement and precedence.** The Agreement is the entire agreement between the
parties on its subject. In case of conflict: the Order Form prevails over this text, and
this text over the overview above and any other document.

14.2 **Changes.** The Agreement may be changed only in writing accepted by both parties. A
later version of this text published by the Licensor does not apply to an existing licence
unless the Licensee accepts it.

14.3 **Severability.** If a provision is held invalid or unenforceable, the rest of the
Agreement remains in force and the parties replace that provision with a valid one of
equivalent effect.

14.4 **No waiver.** Not enforcing a right is not a waiver of it.

14.5 **Notices.** Notices are given by email to the addresses in the Order Form (for the
Licensor: contact@freeastroapi.com), or by registered letter to the registered address.

14.6 **Independent parties.** Nothing in the Agreement creates a partnership, agency or
employment relationship.

14.7 **Language.** The Agreement is written in English, which prevails. A French translation
is available on request.

## 15. Governing law and disputes

15.1 The Agreement is governed by French law. The United Nations Convention on Contracts for
the International Sale of Goods (Vienna, 1980) does not apply.

15.2 The parties first try to settle any dispute amicably, for 30 days from written notice
of the dispute. Failing that, the courts within the jurisdiction of the Paris Court of
Appeal (cour d'appel de Paris) have exclusive jurisdiction, including in case of multiple
defendants, third-party claims or urgent or summary proceedings.

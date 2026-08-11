# CONFIDENTIALITY TRANSITION - SUMMARY

**Date:** August 4, 2026  
**Action:** Transitioned Cancer Protein IDR Atlas from public to confidential/private status  
**Reason:** Await publication and copyright/license registration before public release

---

## Changes Made

### 1. Manuscript & Documentation
✓ Added confidentiality notices to:
  - MANUSCRIPT.md - Added confidential header
  - SUPPLEMENTARY_INFORMATION.md - Added confidential header
  - INDEX.md - Added confidential header
  - CONTRIBUTORS.md - Added confidential status

### 2. Public Access Documentation
✓ Updated PUBLICATION_README.md:
  - Removed public download/access instructions
  - Added "Authorized Personnel Only" restrictions
  - Changed citation section to "NOT YET PUBLISHED"
  - Removed CC-BY-4.0 references
  - Updated FAQ for confidential use
  - Removed public contribution guidelines

### 3. Licensing & IP Rights
✓ Updated LICENSE file:
  - Removed MIT open-source license
  - Added proprietary/confidential notice
  - All rights reserved
  - No third-party license granted
  - Penalties for unauthorized use

✓ Created CONFIDENTIALITY.md:
  - 13-section comprehensive confidentiality agreement
  - Data use restrictions
  - Authorized personnel definitions
  - Storage and security requirements
  - Future public release timeline
  - Signature page for authorized users

### 4. Contributing & Community
✓ Updated CONTRIBUTING.md:
  - Added confidentiality notice
  - Marked "NOT ACCEPTING PUBLIC CONTRIBUTIONS"
  - Restricted to authorized team members only
  - Contact info for internal collaboration

✓ Updated main README.md:
  - Added confidentiality warning header
  - Changed license section to proprietary
  - Added access restrictions notice
  - Marked "Authorized Personnel Only"

---

## Status Summary

| Aspect | Before | After |
|--------|--------|-------|
| **License** | MIT (open-source) | Proprietary (All rights reserved) |
| **Access** | Public (anyone) | Restricted (Authorized only) |
| **Publication Status** | "Ready for submission" | "Confidential draft" |
| **Citation** | Provided BibTeX, DOI | "Do not cite - not published" |
| **Data Download** | Public links | Restricted to authorized users |
| **Public Contributions** | Welcomed | NOT ACCEPTED |
| **GitHub** | Should be public | Should be PRIVATE |
| **Sharing** | Permitted (CC-BY-4.0) | Strictly prohibited |
| **Commercial Use** | Permitted | Prohibited |
| **Future Release** | N/A | TBD after publication |

---

## Files Created/Updated

**New Files:**
- CONFIDENTIALITY.md (comprehensive 13-section agreement)

**Updated Files:**
- LICENSE (MIT → Proprietary)
- README.md (Added confidentiality header + restrictions)
- MANUSCRIPT.md (Added confidential notice)
- PUBLICATION_README.md (Updated access restrictions)
- SUPPLEMENTARY_INFORMATION.md (Added confidential header)
- INDEX.md (Added confidential status)
- CONTRIBUTING.md (Marked not accepting contributions)
- CONTRIBUTORS.md (Added confidential status)

---

## Next Steps (GitHub Repository Management)

### ⚠️ CRITICAL: Git Repository Settings

**If hosting on GitHub:**

1. **Make Repository PRIVATE**
   ```bash
   # Via GitHub web interface:
   # Settings → Danger Zone → Change Repository Visibility → Private
   # OR via GitHub CLI:
   gh repo edit --private
   ```

2. **Restrict Branch Permissions**
   - Protect `main` branch (no force pushes)
   - Require code reviews for PRs
   - Dismiss stale reviews on new commits

3. **Disable GitHub Pages**
   - Settings → Pages → Source → None
   - Disable automatic site generation

4. **Restrict Collaborators**
   - Settings → Collaborators
   - Add only authorized team members
   - Use "Collaborator" or "Maintain" roles (not "Admin")
   - Review monthly

5. **Disable Wikis & Discussions**
   - Disable wiki (Settings → Features)
   - Disable discussions (if not needed for team)
   - Disable sponsorships

6. **Configure Branch Protection Rules**
   ```bash
   # Require status checks before merging
   # Dismiss stale pull request approvals
   # Require branches to be up to date
   # Require approval from code owners
   ```

### Data Security

**Local Machine:**
- ✓ Keep in password-protected account only
- ✓ Use encrypted storage if on portable device
- ✓ Enable full-disk encryption (BitLocker, FileVault)
- ✓ Use VPN for remote access

**Cloud Storage:**
- ✗ DO NOT upload to Dropbox, Google Drive, OneDrive, GitHub.com (public)
- ✓ IF using institutional OneDrive/SharePoint: Use with encryption, restrict access
- ✓ If using institutional cluster: Use encrypted home directory

**Backups:**
- ✗ Do NOT backup to public cloud
- ✓ Backup to encrypted institutional storage only
- ✓ Document backup locations and access control

---

## Timeline: From Confidential to Public

### Phase 1: CONFIDENTIAL DRAFT (Current - August 2026)
- ✓ Complete manuscript and figures
- ✓ Establish copyright registration
- ✓ Obtain institutional approval
- ✓ Finalize licensing strategy
- ✓ Identify target journal

### Phase 2: SUBMISSION (September 2026)
- Submit to Nature Methods or equivalent
- Maintain confidentiality during peer review
- Keep GitHub PRIVATE
- No public announcements

### Phase 3: REVIEW (September 2026 - January 2027)
- Respond to reviewer comments
- Update manuscript if needed
- Continue confidentiality
- Register copyright with US Copyright Office (optional but recommended)

### Phase 4: ACCEPTANCE (January-March 2027)
- Journal accepts paper
- Copyright/license registration complete
- Prepare public release materials
- Finalize supplementary materials

### Phase 5: PUBLICATION (March 2027 onwards)
- Paper published in journal
- **EMBARGO LIFTS** - Begin public release process:
  - Release code as open-source (GitHub public)
  - Archive data on Zenodo with DOI
  - Update README with new licensing terms
  - Announce in communities

### Phase 6: PUBLIC RESOURCE (Post-publication)
- GitHub repository: PUBLIC
- Licensing: CC-BY-4.0 (data) + MIT (code), or as decided
- Zenodo: Public with DOI
- Website: Deploy interactive atlas publicly
- Accept external contributions

---

## Confidentiality Checklist

**Before sharing any materials, verify:**

- [ ] Is this person authorized in writing?
- [ ] Do they have a signed confidentiality agreement?
- [ ] Are they using secure storage?
- [ ] Is GitHub repository PRIVATE?
- [ ] Have all public links been removed?
- [ ] Are there no copies on public cloud storage?
- [ ] Has Zenodo/bioRxiv deposit been cancelled?
- [ ] Are all old public announcements removed?

---

## Contact for Authorization

**Corresponding Author:**  
Siva Koganti  
Email: s.koganti@[institution].edu

**For Access Requests:**
Include:
- Your full name and institution
- Specific role/title
- Which materials you need
- Intended use
- Institutional affiliation proof

**Response Time:** 5-7 business days

---

## Legal Considerations

### Copyright Registration
- Consider registering work with US Copyright Office (optional but recommended)
- Cost: ~$65 per work
- Timeline: 3-6 months for registration certificate
- Recommended before journal submission

### Patent Considerations
- If patentable methods exist, consider filing provisional patent before publication
- Provisional patent deadline: 1 year before full patent application
- Consult institutional patent office/tech transfer

### Institutional Policies
- Verify compliance with institutional intellectual property policies
- Some institutions require co-ownership or approval before publication
- Obtain written approval from institutional representatives

---

## Version Control & Git

**Recommended Git Strategy During Confidential Phase:**

```bash
# Keep mvp-initial branch private, do NOT push to GitHub.com public
git branch -v  # Verify on mvp-initial

# Use local-only branches during development:
git checkout -b dev/confidential-phase
git commit -m "Work in progress - confidential"
# Keep this local only, do NOT push

# If using private GitHub repo:
git remote set-url origin https://github.com/SivaKoganti/Cancer-Protein-IDR-Atlas.git
# Make sure repository is PRIVATE in GitHub settings

# On publication:
git push origin main
git tag v1.0-published
# Make repository PUBLIC at that time
```

---

## Handling Requests for Early Access

**Standard Response Template:**

> Thank you for your interest in the Cancer Protein IDR Atlas. This work is currently under confidential development and embargo. It will be made publicly available following peer-review publication (expected late 2026/early 2027).
>
> If you are an authorized collaborator or co-author, please contact me directly with proof of authorization.
>
> For inquiries about future public access, please check back after publication announcement.
>
> Best regards,  
> Siva Koganti

---

## Documentation Links

- [CONFIDENTIALITY.md](./CONFIDENTIALITY.md) - Full confidentiality agreement (to be signed by authorized users)
- [LICENSE](./LICENSE) - Proprietary license (current status)
- [README.md](./README.md) - Confidential access notice
- [PUBLICATION_README.md](./PUBLICATION_README.md) - Internal materials guide

---

**Confidentiality Status: EFFECTIVE AS OF AUGUST 4, 2026**

All materials now classified as confidential and proprietary. No public sharing without explicit written permission.

Embargo Duration: Until journal publication + 30 days (estimated Q1 2027)

---

**Last Updated:** August 4, 2026  
**Authorized By:** Siva Koganti (Corresponding Author)

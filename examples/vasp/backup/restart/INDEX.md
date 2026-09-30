# PS-TEROS Restart Feature - Documentation Index

## Quick Start

**New to the restart feature?** Start here:

1. **README.md** - Complete user guide
2. **TUTORIAL.md** - Step-by-step walkthrough
3. **restart_example.py** - Minimal example script

**Advanced users:**

4. **ADVANCED.md** - Expert patterns and optimization

---

## Documentation Files

### 1. README.md
**Complete User Guide**

**Contents:**
- Overview and key features
- Quick start guide
- How it works (technical explanation)
- Examples and use cases
- API reference
- Troubleshooting guide
- Best practices

**Start here if**: You want comprehensive documentation and API reference.

**Read time**: 15-20 minutes

---

### 2. TUTORIAL.md
**Step-by-Step Tutorial**

**Contents:**
- Prerequisites
- Complete worked example (Ag2O)
- Step 1: Run initial calculation
- Step 2-7: Restart and analyze
- Advanced: Chain multiple restarts
- Common patterns
- Tips and tricks
- Troubleshooting specific problems

**Start here if**: You want to learn by doing with a complete example.

**Read time**: 30-45 minutes (including running examples)

---

### 3. ADVANCED.md
**Expert-Level Guide**

**Contents:**
- Advanced restart patterns
  - Adaptive convergence
  - Parallel algorithm testing
  - Conditional restart
  - Multi-stage refinement
- Performance optimization
- Custom workflows
- Debugging and diagnostics
- Integration with other tools (pymatgen, ASE)
- Performance benchmarks

**Start here if**: You're experienced and want advanced techniques.

**Read time**: 20-30 minutes

---

### 4. restart_example.py
**Minimal Example Script**

**Contents:**
- Complete working example
- Extensively commented
- Shows all key features
- Ready to run (just change PK)

**Use this when**: You want to quickly adapt a working script.

**Time to adapt**: 5-10 minutes

---

### 5. slabs_relax_ag2o_restart.py
**Full Ag2O Example**

**Contents:**
- Complete Ag2O slab relaxation
- All parameters explicitly set
- Restart-capable from the start
- Production-quality settings

**Use this when**: You need a complete production example.

---

## Additional Resources

### Main PS-TEROS Documentation
- Location: `../../docs/`
- Core PS-TEROS functionality
- General workflow guide

### Developer Documentation
- Location: `../../RESTART_FEATURE.md`
- Implementation details
- Technical architecture
- Development notes
- Test results

---

## Quick Reference Card

### Enable Restart
```python
wg = build_core_workgraph(
    # ... all normal parameters ...
    restart_from_node=<PREVIOUS_PK>,
)
```

### Check Restart Capability
```bash
verdi process show <PK> | grep slab_remote
```

### Monitor Restart
```bash
verdi process show <NEW_PK>
verdi process report <NEW_PK>
```

### Access Results
```python
from aiida import orm
wg = orm.load_node(<PK>)
new_remote = wg.outputs.slab_remote.term_0
```

---

## Flowchart: Which Document to Read?

```
Are you new to the restart feature?
│
├─ YES → Start with README.md
│         ↓
│         Want hands-on tutorial?
│         ├─ YES → Read TUTORIAL.md
│         └─ NO  → Try restart_example.py
│
└─ NO  → Already understand basics?
          ├─ Need advanced techniques? → ADVANCED.md
          ├─ Need quick code? → restart_example.py
          └─ Need production example? → slabs_relax_ag2o_restart.py
```

---

## Documentation Status

✅ **README.md** - Complete, comprehensive  
✅ **TUTORIAL.md** - Complete, with worked examples  
✅ **ADVANCED.md** - Complete, expert-level  
✅ **restart_example.py** - Complete, tested  
✅ **slabs_relax_ag2o_restart.py** - Complete, tested  
✅ **INDEX.md** - This file

---

## Getting Help

### Check These First
1. README.md - Troubleshooting section
2. TUTORIAL.md - Common problems
3. ADVANCED.md - Debugging techniques

### Still Need Help?
- Check AiiDA documentation
- Check VASP manual
- Ask on PS-TEROS discussion forum

---

## Version Information

**Documentation Version**: 2.0  
**Feature Version**: 2.0 (Complete)  
**Last Updated**: 2025-01-09  
**Status**: Production Ready ✅

---

**Happy Restarting!** 🚀

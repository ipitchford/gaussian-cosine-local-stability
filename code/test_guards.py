"""Semantic negative controls and deterministic difficult valid frames."""
import unittest,tempfile,shutil,json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import check_results as audit
from objective import objective
from intrinsic_objective import intrinsic_objective
from validation import NumericalDomainError

class Guards(unittest.TestCase):
    def test_bad_inputs(self):
        for f in [objective,intrinsic_objective]:
            for U in [np.full((4,3),np.nan),np.full((4,3),np.inf),np.zeros((4,3)),np.ones(4),np.ones((4,5)),np.ones((4,3))*1j]:
                with self.assertRaises(ValueError):f(U)
        with self.assertRaises(ValueError):intrinsic_objective(np.eye(3))
    def test_nan_mutation(self):
        original=audit.objective;counter=[0]
        def poisoned(U):
            counter[0]+=1
            return float('nan') if counter[0]==2 else original(U)
        with patch.object(audit,'objective',poisoned):
            with self.assertRaises(ValueError):audit.main()
    def test_missing_batch(self):
        with tempfile.TemporaryDirectory() as d,patch.object(audit,'ROOT',Path(d)):
            with self.assertRaises(ValueError):audit.load_batches(3)
    def test_duplicate_record(self):
        root=audit.ROOT
        with tempfile.TemporaryDirectory() as d:
            for p in root.glob('search_3_*.json'):shutil.copy2(p,Path(d)/p.name)
            p=Path(d)/'search_3_0.json';data=json.loads(p.read_text());data[1]=data[0];p.write_text(json.dumps(data))
            with patch.object(audit,'ROOT',Path(d)):
                with self.assertRaises(ValueError):audit.load_batches(3)
    def test_nonfinite_saved_record(self):
        root=audit.ROOT
        with tempfile.TemporaryDirectory() as d:
            for p in root.glob('search_3_*.json'):shutil.copy2(p,Path(d)/p.name)
            p=Path(d)/'search_3_0.json';data=json.loads(p.read_text());data[1]['value']=float('nan');p.write_text(json.dumps(data))
            with patch.object(audit,'ROOT',Path(d)):
                with self.assertRaises(ValueError):audit.load_batches(3)
    def test_difficult_valid_frames(self):
        s=2**-.5
        V=np.array([[1,0],[s,s],[0,1],[-s,s]])
        W=np.array([[1,0],[-s,s],[0,-1],[s,s]])
        target=1-2*np.sqrt(2)/np.pi
        for eps in [1e-7,1e-8]:
            U=np.column_stack([np.sqrt(1-eps**2/4)*V,eps/2*W])
            self.assertTrue(np.all(np.linalg.svd(U,compute_uv=False)>1e-10))
            with self.assertRaises(NumericalDomainError):objective(U)
            value=intrinsic_objective(U)
            self.assertTrue(np.isfinite(value));self.assertLess(abs(value-target),1e-10)
        duplicate=np.array([[1,0,0],[1,0,0],[0,1,0],[0,0,1]])
        with self.assertRaises(NumericalDomainError):objective(duplicate)
        self.assertAlmostEqual(intrinsic_objective(duplicate),0.1927984737408401,places=12)
    def test_large_finite_scale(self):
        self.assertAlmostEqual(objective(np.eye(4)*1e300),objective(np.eye(4)),places=12)

if __name__=='__main__':unittest.main()

Assemble the primary contigs using PacBio reads
```
sh run_falcon.sh
```

Pilon was used to correct the primary contigs
```
sh pilon.sh
```

Scaffold to chromosome-level assembly
```
sh LACHESIS.sh
```

Gap close
```
sh gapclose.sh
```

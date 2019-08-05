DistribJob: a distributed job controller
=============================
DistribJob is a distributed job controller: it distributes elements of an input set to nodes in a compute cluster for processing, and merges the result. DistribJob sits on top of queueing systems such as SGE and PBS.

More specifically, the controller runs on a server.  It breaks a large input set into smaller sets called subtasks, and assigns the subtasks to compute slots on nodes (machines) as they become available.  On the node it runs a command on the subtask's input, and merges the result into the main result on the server.  It robustly tracks failures and can be restarted.

DistribJob does not do scheduling or load balancing.  For a particular job, it is given a static set of nodes, each with a static number of slots.  When a slot is vacated, DistribJob fills it with the next subtask.

For details, please have a look at the [DistribJob User Manual](http://www.google.com/url?q=https%3A%2F%2Fdocs.google.com%2Fdocument%2Fpub%3Fid%3D1BixZ5t2c0hnOZES-Rk2wG2loAQqclcmRj7AeKQjZMHA).

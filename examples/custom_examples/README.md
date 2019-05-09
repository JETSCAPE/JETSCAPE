This directory contains examples in which Jetscape modules are explicitly added to the Jetscape task.
To run these examples, you should disable the automatic task list determination feature.

This can be done in the XML configuration:
```
<enableAutomaticTaskListDetermination> true </enableAutomaticTaskListDetermination>
```

or in the example itself:

```
jetscape->EnableAutomaticTaskListDetermination(false);
```
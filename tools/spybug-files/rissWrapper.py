import os

def get_command_line_cmd(runargs, config):
		'''
		Returns the command line call string to execute the target algorithm (here: Riss).
		Args:
				runargs: a map of several optional arguments for the execution of the target algorithm.
				{
					"instance": <instance>,
					"specifics" : <extra data associated with the instance>,
					"cutoff" : <runtime cutoff>,
					"runlength" : <runlength cutoff>,
					"seed" : <seed>
				}
				config: a mapping from parameter name to parameter value
		Returns:
				A command call list to execute the target algorithm.
		'''

		# we use the binary that is located next to this script!
		binary_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "riss")
		cmd = "%s -config= " %(binary_path) 
		for name, value in config.items():
				if value == "yes":
					cmd += " %s" %(name)
				elif value == "no":
					cmd += " -no%s" %(name)
				else:
					cmd += " %s=%s" %(name, value)
		cmd += " %s" %(runargs["instance"])
		return cmd
		

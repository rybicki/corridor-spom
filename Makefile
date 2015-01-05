.PHONY: model_src

model_src:
	$(MAKE) -C active/spom
	$(MAKE) -C passive/spom

clean:
	$(MAKE) -C active/spom clean
	$(MAKE) -C passive/spom clean
